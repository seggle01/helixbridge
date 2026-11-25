"""Generate a target-only MSA and save it as target.a3m."""
import argparse
import json
import tarfile
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


def read_json(request):
    for attempt in range(3):
        try:
            with urlopen(request, timeout=120) as response:
                return json.load(response)
        except TimeoutError:
            if attempt == 2:
                raise RuntimeError(
                    "MMseqs MSA service is unreachable at "
                    "https://api.colabfold.com after 3 attempts"
                )
            print("Connection timed out; retrying...", flush=True)
            time.sleep(10)


def read_first_fasta_sequence(path):
    records = []
    name = None
    residues = []

    for line in path.read_text().splitlines() + [">END"]:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(residues)))
            name = line[1:]
            residues = []
        else:
            residues.append(line)

    if not records or not records[0][1]:
        raise ValueError("FASTA does not contain a non-empty first record")
    return records[0][1]


parser = argparse.ArgumentParser(description="Generate a target-only MSA from a FASTA file")
parser.add_argument("--fasta", required=True, help="Input FASTA file")
parser.add_argument("--output", required=True, help="Output directory")
args = parser.parse_args()

sequence = read_first_fasta_sequence(Path(args.fasta))
output = Path(args.output)
output.mkdir(parents=True, exist_ok=True)

server = "https://api.colabfold.com"
headers = {"User-Agent": "helixbridge/1.0"}
payload = urlencode({"q": ">query\n" + sequence + "\n", "mode": "env"}).encode()
reply = read_json(Request(server + "/ticket/msa", data=payload, headers=headers))
ticket_id = reply.get("id")
if not ticket_id:
    raise RuntimeError("MSA service did not accept the request: " + json.dumps(reply))

deadline = time.monotonic() + 600
while True:
    reply = read_json(Request(server + "/ticket/" + ticket_id, headers=headers))
    status = reply.get("status")
    print("MSA state:", status, flush=True)
    if status == "COMPLETE":
        break
    if status in ("ERROR", "MAINTENANCE", "UNKNOWN", "RATELIMIT"):
        raise RuntimeError("MSA service status: " + status)
    if time.monotonic() > deadline:
        raise TimeoutError("MSA search timed out")
    time.sleep(10)

archive = output / "result.tar.gz"
with urlopen(Request(server + "/result/download/" + ticket_id, headers=headers), timeout=90) as response:
    archive.write_bytes(response.read())

entries = []
with tarfile.open(archive) as tar:
    for member in tar.getmembers():
        if not member.name.endswith(".a3m") or not member.isfile():
            continue
        if member.size > 100_000_000:
            raise ValueError("Unexpectedly large MSA result")
        raw = tar.extractfile(member).read().decode().replace("\x00", "")
        name = None
        seq = []
        for line in raw.splitlines() + [">END"]:
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                if name is not None:
                    entries.append((name, "".join(seq)))
                name = line[1:]
                seq = []
            else:
                seq.append(line.strip())

seen = set()
kept = []
for name, seq in entries:
    if not seq:
        continue
    aligned = "".join(character for character in seq if not character.islower())
    if len(aligned) != len(sequence):
        raise ValueError("MSA alignment width mismatch")
    if seq not in seen:
        seen.add(seq)
        kept.append((name, seq))

if not kept or kept[0][1] != sequence:
    raise ValueError("MSA query does not match the FASTA target")

kept = kept[:1024]
if len(kept) < 2:
    raise ValueError("No homologs returned")

path = output / "target.a3m"
path.write_text("".join(">" + name + "\n" + seq + "\n" for name, seq in kept))
archive.unlink()
print("Saved target-only MSA: " + str(len(kept)) + " rows")
