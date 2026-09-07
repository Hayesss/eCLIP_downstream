"""Download and verify the fixed GENCODE v38 / GRCh38 reference resources."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
BASE_URL = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/'
NAME = 'GRCh38.primary_assembly.genome.fa.gz'
EXPECTED_SIZE = 844691642
EXPECTED_MD5 = '71b78189dc2a7b7e7af802949f168aa7'


def digest(path, algorithm='sha256'):
    value = hashlib.new(algorithm)
    with path.open('rb') as stream:
        while block := stream.read(8 * 1024 * 1024):
            value.update(block)
    return value.hexdigest()


def download_file(url, destination):
    temporary = destination.with_name(destination.name + '.partial')
    with urllib.request.urlopen(url, timeout=60) as source, temporary.open('wb') as out:
        shutil.copyfileobj(source, out, 1024 * 1024)
    temporary.replace(destination)


def download_ranges(url, destination, size, workers=8, block_size=8*1024*1024):
    """Resume fixed file-byte chunks; reject servers that ignore Range headers."""
    if workers < 1 or block_size < 1 or size < 1:
        raise ValueError('workers, block_size and size must be positive')
    parts = destination.parent / (destination.name + '.parts')
    parts.mkdir(exist_ok=True)

    def fetch(start):
        end = min(size, start+block_size)-1
        part = parts / f'{start:012d}.part'
        if part.exists() and part.stat().st_size == end-start+1:
            return part
        for attempt in range(3):
            try:
                request = urllib.request.Request(url, headers={'Range': f'bytes={start}-{end}'})
                temporary = part.with_suffix('.partial')
                with urllib.request.urlopen(request, timeout=60) as source:
                    expected_range = f'bytes {start}-{end}/{size}'
                    if source.status != 206 or source.headers.get('Content-Range') != expected_range:
                        raise ValueError('Server did not return the requested byte range')
                    with temporary.open('wb') as out:
                        shutil.copyfileobj(source, out, 1024*1024)
                if temporary.stat().st_size != end-start+1:
                    raise ValueError('Incomplete reference download chunk')
                temporary.replace(part)
                return part
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(attempt+1)

    starts = list(range(0, size, block_size))
    complete = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for job in as_completed([pool.submit(fetch, start) for start in starts]):
            complete += job.result().stat().st_size
            print(f'Downloaded {complete/1e6:.1f}/{size/1e6:.1f} MB', flush=True)
    temporary = destination.with_name(destination.name + '.assembled.partial')
    with temporary.open('wb') as out:
        for start in starts:
            with (parts / f'{start:012d}.part').open('rb') as source:
                shutil.copyfileobj(source, out, 8*1024*1024)
    temporary.replace(destination)
    return parts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'references'/'gencode_v38')
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--metadata-only', action='store_true', help='Download only official checksums and gene/transcript ID lists')
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers must be positive')
    base = args.output_dir.resolve()
    base.mkdir(parents=True, exist_ok=True)
    md5_file = base/'MD5SUMS'
    download_file(BASE_URL+'MD5SUMS', md5_file)
    checks = {line.split()[1].lstrip('./'): line.split()[0] for line in md5_file.read_text().splitlines() if len(line.split()) == 2}
    if checks[NAME] != EXPECTED_MD5:
        raise ValueError('Official reference checksum differs from the pinned GENCODE v38 reference')
    metadata = {}
    for kind in ['Gene', 'Transcript']:
        name = f'gencode.v38.metadata.{kind}_source.gz'
        path = base/name
        if not path.exists() or digest(path, 'md5') != checks[name]:
            download_file(BASE_URL+name, path)
        if digest(path, 'md5') != checks[name]:
            raise ValueError(f'MD5 mismatch: {name}')
        metadata[name] = {'url': BASE_URL+name, 'md5': checks[name], 'sha256': digest(path)}
    (base/'metadata_manifest.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    if args.metadata_only:
        print('GENCODE v38 ID metadata verified')
        return
    archive = base/NAME
    parts = None
    if not archive.exists():
        parts = download_ranges(BASE_URL+NAME, archive, EXPECTED_SIZE, args.workers)
    if archive.stat().st_size != EXPECTED_SIZE or digest(archive, 'md5') != EXPECTED_MD5:
        raise ValueError(f'Reference checksum mismatch. Remove the invalid cache file and its .parts directory, then retry: {archive}')
    print('Official reference MD5 verified', flush=True)
    fasta = base/NAME[:-3]
    manifest_path = base/'reference_manifest.json'
    previous = json.loads(manifest_path.read_text(encoding='utf-8')) if manifest_path.exists() else {}
    trusted_fasta = fasta.exists() and previous.get('sha256_fasta') and digest(fasta) == previous['sha256_fasta'] and previous.get('md5') == EXPECTED_MD5
    if not trusted_fasta:
        print('Decompressing verified reference', flush=True)
        temporary = fasta.with_name(fasta.name+'.partial')
        with gzip.open(archive, 'rb') as source, temporary.open('wb') as out:
            shutil.copyfileobj(source, out, 8*1024*1024)
        temporary.replace(fasta)
        # An index associated with a replaced FASTA must be regenerated.
        fasta.with_name(fasta.name+'.fai').unlink(missing_ok=True)
    report = {'source_url': BASE_URL+NAME, 'gencode_release': 38, 'assembly': 'GRCh38 primary assembly', 'md5': EXPECTED_MD5, 'sha256_gz': digest(archive), 'sha256_fasta': digest(fasta), 'fasta_path': str(fasta), 'compressed_bytes': archive.stat().st_size, 'fasta_bytes': fasta.stat().st_size, 'metadata': metadata}
    manifest_path.write_text(json.dumps(report, indent=2), encoding='utf-8')
    if parts is not None:
        for part in parts.iterdir():
            if part.is_file() and part.suffix in ['.part', '.partial']:
                part.unlink()
        parts.rmdir()
    print(f'Reference ready: {fasta}', flush=True)


if __name__ == '__main__':
    main()
