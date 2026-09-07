"""Deterministic tests using synthetic DNA; no private data or network access."""
import io
from pathlib import Path
import random
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from extract_sequences import Genome, interval, reverse_complement
from download_reference import download_ranges


class SequenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        rng = random.Random(38)
        self.sequences = {name: ''.join(rng.choice('ACGT') for _ in range(size)) for name,size in [('chrTest',1800),('chrM',16569)]}
        self.path = Path(self.temp.name)/'synthetic.fa'
        self.path.write_text(''.join('>'+name+' description\n'+'\n'.join(seq[i:i+60] for i in range(0,len(seq),60))+'\n' for name,seq in self.sequences.items()), encoding='ascii')
        self.genome = Genome(self.path)

    def tearDown(self):
        self.genome.data.close()
        self.genome.stream.close()
        self.temp.cleanup()

    def test_index_against_known_sequences(self):
        rng = random.Random(500)
        for name, seq in self.sequences.items():
            self.assertEqual(self.genome.index[name][0], len(seq))
            for _ in range(200):
                start = rng.randrange(len(seq))
                end = rng.randrange(start, len(seq)+1)
                self.assertEqual(self.genome.fetch(name,start,end), seq[start:end])

    def test_strand_windows_and_chromosome_boundaries(self):
        for name,seq in self.sequences.items():
            for strand in ['+','-']:
                for start,end in [(0,20),(790,851),(len(seq)-20,len(seq))]:
                    record = {'chr':name,'start':start,'end':end,'strand':strand}
                    for target in [500,1000]:
                        _,_,segments,offset,boundary = interval(record,target,self.genome)
                        result = ''.join(seq[a:b] for a,b in segments)
                        peak = seq[start:end]
                        if strand == '-':
                            result,peak = reverse_complement(result),reverse_complement(peak)
                        self.assertEqual(len(result),target)
                        self.assertEqual(result[offset:offset+end-start],peak)
                        if start == 0 or end == len(seq):
                            expected = 'circular_wrap_chrM' if name == 'chrM' else 'shifted_inside_chromosome'
                            self.assertEqual(boundary,expected)

    def test_odd_flank_base_goes_downstream(self):
        for strand in ['+','-']:
            record = {'chr':'chrTest','start':790,'end':851,'strand':strand}
            _,_,_,offset,_ = interval(record,500,self.genome)
            self.assertEqual(offset,219)
            self.assertEqual(500-offset-61,220)

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement('AACG'),'CGTT')
        self.assertEqual(reverse_complement('ACGTN'),'NACGT')

    def test_reject_peak_larger_than_window(self):
        with self.assertRaises(AssertionError):
            interval({'chr':'chrTest','start':100,'end':700,'strand':'+'},500,self.genome)


class RangeDownloadTests(unittest.TestCase):
    def test_ranges_and_cached_resume(self):
        payload = bytes(range(251))*3
        calls = []
        def response(request, timeout):
            begin,end = map(int,request.get_header('Range').split('=')[1].split('-'))
            calls.append((begin,end))
            result = io.BytesIO(payload[begin:end+1])
            result.status = 206
            result.headers = {'Content-Range':f'bytes {begin}-{end}/{len(payload)}'}
            return result
        with tempfile.TemporaryDirectory() as folder, patch('download_reference.urllib.request.urlopen',side_effect=response):
            path = Path(folder)/'reference.gz'
            download_ranges('https://example.invalid/reference',path,len(payload),2,256)
            self.assertEqual(path.read_bytes(),payload)
            self.assertEqual(sorted(calls),[(0,255),(256,511),(512,752)])
            calls.clear()
            download_ranges('https://example.invalid/reference',path,len(payload),2,256)
            self.assertEqual(path.read_bytes(),payload)
            self.assertEqual(calls,[])

    def test_server_ignoring_range_is_rejected(self):
        def response(*args, **kwargs):
            result = io.BytesIO(b'abc')
            result.status = 200
            result.headers = {}
            return result
        with tempfile.TemporaryDirectory() as folder, patch('download_reference.urllib.request.urlopen',side_effect=response), patch('download_reference.time.sleep'):
            path = Path(folder)/'reference.gz'
            with self.assertRaisesRegex(ValueError,'byte range'):
                download_ranges('https://example.invalid/reference',path,3,1,3)
            self.assertFalse(path.exists())


if __name__ == '__main__':
    unittest.main()
