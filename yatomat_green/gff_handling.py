from typing import List, Dict
from dataclasses import dataclass

@dataclass
class GenomeFeature:
    """Class for storing genomic feature information"""
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str = "."
    strand: str = "+"
    phase: str = "."
    attributes: Dict[str, str] = None

# Remove GFFWriter class from here


class GFFWriter:
    """Utility class for writing GFF3 files"""
    
    @staticmethod
    def format_attributes(attrs: Dict[str, str]) -> str:
        if not attrs:
            return "."
        return ";".join(f"{k}={v}" for k, v in attrs.items())
    
    def write(self, features: List[GenomeFeature], output_path: str):
        with open(output_path, 'w') as f:
            f.write("##gff-version 3\n")
            for feature in features:
                attrs = self.format_attributes(feature.attributes or {})
                f.write(f"{feature.seqid}\t{feature.source}\t{feature.type}\t"
                        f"{feature.start}\t{feature.end}\t{feature.score}\t"
                        f"{feature.strand}\t{feature.phase}\t{attrs}\n")
