import json
import hashlib
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[1]
class DeliveryTests(unittest.TestCase):
    def test_scientific_claims_and_images(self):
        n=json.loads((ROOT/'三体问题数值求解.ipynb').read_text())
        text='\n'.join(''.join(c['source']) for c in n['cells'])
        for bad in ['数值解正确 ⟺','微扰应只带来有界偏差','尚处指数增长的早期','反复弹甩，轨迹混乱但被束缚','<img src="<img','效率提升 3.5 倍','没有相位漂移导致的涂抹']:
            self.assertNotIn(bad,text)
        self.assertIn('delta_v2 = 1e-6',text)
        self.assertIn('figure8_validation',text)
        self.assertIn('basic_validation',text)
    def test_reference_provenance(self):
        meta=json.loads((ROOT/'data/threebody_reference/metadata.json').read_text())
        self.assertEqual(meta['method'],'DOP853')
        for name, digest in meta['files'].items():
            self.assertEqual(hashlib.sha256((ROOT/'data/threebody_reference'/name).read_bytes()).hexdigest(),digest)
        self.assertEqual(set(meta['cases']),{'single','two','three','figure8'})
    def test_validation_source_matches_notebook(self):
        n=json.loads((ROOT/'三体问题数值求解.ipynb').read_text())
        cells=[''.join(c['source']) for c in n['cells']]
        self.assertIn((ROOT/'src/threebody/validation.jl').read_text(),cells)
if __name__=='__main__':unittest.main()
