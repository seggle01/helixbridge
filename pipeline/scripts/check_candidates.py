"""Contest-readiness check for candidate binders: novelty and duplicate detection.

Applies the two competition rules to an RFD4-Proteina candidate FASTA:

NOVELTY   a candidate is too similar to a known binder when sequence identity
          is > 70% over an alignment covering >= 80% of the full length. Both
          the candidate and the known binder must meet the coverage gate: a
          short fragment would otherwise "cover" 100% of itself while aligning
          to a small corner of a much larger protein.

DUPLICATE a later sequence duplicates an earlier one when it is 100% identical
          over >= 90% of the earlier sequence. For antibodies only the CDR
          regions are compared; minibinders are compared whole.

Identity is measured over the alignment span. Coverage is measured on a
full-length (global) alignment, because the reportable number is identity over
the whole molecule, not over the best local island. A local alignment that
truncates to a shared motif reports high identity over a handful of residues and
would score a ten-residue overlap the same as a 100% identical duplicate.

The duplicate rule gates on coverage of the EARLIER sequence only, matching the
wording of the rule, and requires strictly 100% identity over that overlap.

Usage:
  python check_candidates.py --candidates candidates.fasta --novelty
  python check_candidates.py --candidates candidates.fasta --dedup
  python check_candidates.py --candidates candidates.fasta --novelty --dedup
"""
import argparse
import csv
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

IDENTITY_LIMIT = 0.70
NOVELTY_COVERAGE = 0.80
DUP_COVERAGE = 0.90
POLY_X_WINDOW = 8
POLY_X_FRACTION = 0.5


def build_aligner(mode):
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


GLOBAL = build_aligner('global')
LOCAL = build_aligner('local')


def identity_over_span(aligner, query, target):
    """Identity across the alignment span, with both coverages of that span."""
    try:
        alignment = aligner.align(query, target)[0]
    except IndexError:
        return 0.0, 0.0, 0.0
    identical = 0
    for (qs, qe), (ts, te) in zip(*alignment.aligned):
        for q, t in zip(query[qs:qe], target[ts:te]):
            identical += q == t
    q_start, q_end = alignment.aligned[0][0][0], alignment.aligned[0][-1][1]
    t_start, t_end = alignment.aligned[1][0][0], alignment.aligned[1][-1][1]
    span = max(q_end - q_start, t_end - t_start)
    if span == 0:
        return 0.0, 0.0, 0.0
    return identical / span, (q_end - q_start) / len(query), (t_end - t_start) / len(target)


def low_complexity(sequence):
    """Longest homopolymer run and the largest single-residue fraction."""
    if not sequence:
        return 0, 0.0
    longest = run = 1
    for previous, current in zip(sequence, sequence[1:]):
        run = run + 1 if current == previous else 1
        longest = max(longest, run)
    fraction = max((sequence.count(aa) for aa in set(sequence)), default=0) / len(sequence)
    return longest, fraction


def flag_low_complexity(sequence):
    longest, fraction = low_complexity(sequence)
    reasons = []
    if longest > POLY_X_WINDOW:
        reasons.append(f'homopolymer_run={longest}')
    if fraction > POLY_X_FRACTION:
        reasons.append(f'single_residue_fraction={fraction:.2f}')
    return ';'.join(reasons)


def read_fasta(path):
    records = []
    for record in SeqIO.parse(str(path), 'fasta'):
        residues = str(record.seq).upper()
        if not residues:
            raise SystemExit(f'Empty sequence: {record.id}')
        records.append((record.id.replace(',', '_'), residues))
    if not records:
        raise SystemExit(f'No sequences in {path}')
    return records


def write_csv(path, header, rows):
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def passes_novelty_gate(identity, cand_coverage, known_coverage):
    """Both molecules must be substantially aligned; identity must exceed the limit."""
    return (identity > IDENTITY_LIMIT
            and cand_coverage >= NOVELTY_COVERAGE
            and known_coverage >= NOVELTY_COVERAGE)


def novelty(rows, known):
    header = ['candidate', 'length', 'closest_known', 'identity', 'cand_coverage',
              'known_coverage', 'alignment', 'verdict']
    report = []
    for cid, sequence in rows:
        best = None
        for known_id, known_sequence in known:
            g_id, g_qcov, g_tcov = identity_over_span(GLOBAL, sequence, known_sequence)
            flagged = passes_novelty_gate(g_id, g_qcov, g_tcov)
            # rank flagged pairs first, then by similarity to the aligned region
            key = (flagged, g_id * min(g_qcov, g_tcov))
            if best is None or key > best[0]:
                best = (key, flagged, known_id, g_id, g_qcov, g_tcov)
        _, flagged, known_id, g_id, g_qcov, g_tcov = best
        report.append([cid, len(sequence), known_id, round(g_id, 3),
                       round(g_qcov, 3), round(g_tcov, 3), 'global',
                       'TOO_SIMILAR' if flagged else 'OK'])
    return header, report


def novelty_with_local_crosscheck(rows, known):
    """Count pairs where a local alignment would flag but the global one does not."""
    disagreements = 0
    for cid, sequence in rows:
        for known_id, known_sequence in known:
            g_id, g_qcov, g_tcov = identity_over_span(GLOBAL, sequence, known_sequence)
            l_id, l_qcov, l_tcov = identity_over_span(LOCAL, sequence, known_sequence)
            if passes_novelty_gate(l_id, l_qcov, l_tcov) and not passes_novelty_gate(
                    g_id, g_qcov, g_tcov):
                disagreements += 1
    return disagreements


def dedup(rows):
    header = ['earlier', 'later', 'identity', 'coverage_of_earlier',
              'coverage_of_later', 'span_columns', 'verdict']
    report = []
    for i, (earlier_id, earlier_sequence) in enumerate(rows):
        for later_id, later_sequence in rows[i + 1:]:
            # dup requires a single gapless exact match, so the local alignment
            # is the right tool: a global one spreads a truncated copy across
            # the full length with gaps and reports 96% identity.
            ident, later_cov, earlier_cov = identity_over_span(
                LOCAL, later_sequence, earlier_sequence)
            span = max(round(earlier_cov * len(earlier_sequence)),
                       round(later_cov * len(later_sequence)))
            duplicate = ident >= 0.999 and earlier_cov >= DUP_COVERAGE
            if duplicate:
                report.append([earlier_id, later_id, round(ident, 3),
                               round(earlier_cov, 3), round(later_cov, 3), span,
                               'DUPLICATE'])
    return header, report


def summarize_duplicates(rows, pairs):
    suppresses = {}
    for earlier_id, later_id, _, _, _, _, _ in pairs:
        suppresses.setdefault(later_id, []).append(earlier_id)
    header = ['candidate', 'length', 'low_complexity', 'suppressed_by', 'submit']
    report = []
    for cid, sequence in rows:
        reasons = flag_low_complexity(sequence)
        blockers = suppresses.get(cid, [])
        if blockers:
            submit = 'DROP_DUPLICATE'
        elif reasons:
            submit = 'REVIEW_LOW_COMPLEXITY'
        else:
            submit = 'KEEP'
        report.append([cid, len(sequence), reasons or 'none', ';'.join(blockers) or 'none', submit])
    return header, report


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--candidates', type=Path, required=True)
    parser.add_argument('--known', type=Path, default=None,
                        help='Known-binder FASTA. Default: the IL7Ra corpus beside this script.')
    parser.add_argument('--novelty', action='store_true', help='Run the identity/coverage novelty screen.')
    parser.add_argument('--dedup', action='store_true', help='Run the exact-duplicate scan.')
    parser.add_argument('--out-dir', type=Path, default=Path('.'))
    args = parser.parse_args()
    if not (args.novelty or args.dedup):
        parser.error('choose --novelty and/or --dedup')

    known_path = args.known or (Path(__file__).resolve().parent.parent
                                / 'data/il7ra/il7ra_known_binders/known_il7ra_binders.fasta')
    rows = read_fasta(args.candidates)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.dedup:
        header, pairs = dedup(rows)
        write_csv(args.out_dir / 'dedup_report.csv', header, pairs)
        header, summary = summarize_duplicates(rows, pairs)
        write_csv(args.out_dir / 'duplicate_summary.csv', header, summary)
        kept = sum(1 for row in summary if row[-1] == 'KEEP')
        review = sum(1 for row in summary if row[-1] == 'REVIEW_LOW_COMPLEXITY')
        print(f'candidates        : {len(rows)}')
        print(f'duplicate pairs   : {len(pairs)}')
        print(f'distinct to submit: {len(rows) - len(pairs)}')
        print(f'  ready (KEEP)    : {kept}')
        print(f'  needs review    : {review}')
        for pair in pairs[:10]:
            print(f'  duplicate: {pair[1]} == {pair[0]} '
                  f'(coverage {pair[3]}, identity {pair[2]})')

    if args.novelty:
        known = read_fasta(known_path)
        header, report = novelty(rows, known)
        write_csv(args.out_dir / 'novelty_report.csv', header, report)
        flagged = sum(1 for row in report if row[-1] == 'TOO_SIMILAR')
        identities = sorted(row[3] for row in report)
        margin = IDENTITY_LIMIT - identities[-1]
        print(f'known binders     : {len(known)}')
        print(f'TOO_SIMILAR       : {flagged}')
        print(f'highest identity  : {identities[-1]:.3f} (limit {IDENTITY_LIMIT:.2f}, '
              f'margin {margin:+.3f})')
        print(f'median identity   : {identities[len(identities) // 2]:.3f}')
        print(f'local/global gate disagreements: '
              f'{novelty_with_local_crosscheck(rows, known)}')
    print(f'reports written to: {args.out_dir.resolve()}')


if __name__ == '__main__':
    main()
