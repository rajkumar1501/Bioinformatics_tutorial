import sys
from colored import fg, attr

class Rsid:
    """genetic rsid detector class"""

    def __init__(self, rsid: str, alleles: tuple | str, association: str) -> None:
        """constructor"""
        self.rsid = rsid
        self.alleles = alleles if isinstance(alleles, tuple) else (alleles,)
        self.association = association

    def detect(self, ref: str, alt: str) -> bool:
        """detect this rsid with the given reference and alternate alleles from VCF"""
        print(f"Checking RSID: {self.rsid} | Expected alleles: {self.alleles} | Found alleles: {ref};{alt}")
        detected = any(f"{ref};{alt}" in match for match in self.allele_match_set)
        if detected:
            print(f"Detected RSID {self.rsid} with alleles {ref};{alt}")
        else:
            print(f"RSID {self.rsid} found, but alleles {ref};{alt} do not match expected {self.alleles}")
        return detected

    def __str__(self) -> str:
        """string representation"""
        return f"<Rsid[{self.rsid}]({self.association})>"

    def __repr__(self) -> str:
        """repl representation"""
        return self.__str__()

    @property
    def single_allele(self) -> bool:
        """if this is a single allele match"""
        return len(self.alleles) == 1

    @property
    def allele_match_set(self) -> tuple:
        """form a set of different formats of alleles"""
        def _set(allele: str) -> tuple:
            """form a set of allele formats"""
            if ";" in allele:
                a0, a1 = allele.upper().split(";")
                return (f"{a0};{a1}", f"{a0}|{a1}", f"{a0}{a1}")
            return (allele,)

        tuples = ()
        for allele in self.alleles:
            tuples += _set(allele)
        return tuples

    @property
    def url(self) -> str:
        """return the url to snpedia"""
        return f"https://www.snpedia.com/index.php/{self.rsid}"

    @property
    def info(self) -> str:
        """return an informational string about this rsid when detected"""
        allele = (
            f"with the genotype of ({self.alleles[0]}) " if self.single_allele else ""
        )
        return (
            f"You have a trait which is associated with {trait_colors[self.association]}{self.association}{attr(0)} | "
            f"{self.rsid} {allele}{self.url}"
        )

# Traits mapping to their respective color codes
trait_colors = {
    "Alzheimer's Disease": fg('light_yellow'),
    "Autism": fg('green'),
    "Bipolar Disorder": fg('yellow'),
    "Immunity": fg('red'),
    "Intelligence": fg('cyan'),
    "Longevity": fg('blue'),
    "Metabolism": fg('dark_turquoise'),
    "Muscular Performance": fg('light_green'),
    "OCD": fg('light_magenta'),
    "Schizophrenia": fg('magenta'),
    "Eyes": fg('light_blue'),
    "Hair": fg('light_red'),
}

# Your updated RSID list
_RSIDS = [
    ("rs1006737", "A;A", "Bipolar Disorder"),
    ("rs4027132", "A;A", "Bipolar Disorder"),
    ("rs7570682", "A;A", "Bipolar Disorder"),
    ("rs1375144", "C;C", "Bipolar Disorder"),
    ("rs683395", "C;T", "Bipolar Disorder"),
    ("rs2609653", "C;T", "Bipolar Disorder"),
    ("rs10982256", "C;C", "Bipolar Disorder"),
    ("rs11622475", "C;C", "Bipolar Disorder"),
    ("rs1344484", "T;T", "Bipolar Disorder"),
    ("rs2953145", "C;G", "Bipolar Disorder"),
    ("rs420259", "T;T", "Bipolar Disorder"),
    ("rs4276227", "C;C", "Bipolar Disorder"),
    ("rs4027132", "A;G", "Bipolar Disorder"),
    ("rs2609653", "C;C", "Bipolar Disorder"),
    ("rs2953145", "G;G", "Bipolar Disorder"),
    
    ("rs429358", "C;C", "Alzheimer's Disease"),
    ("rs145999145", "A;A", "Alzheimer's Disease"),
    ("rs908832", "A;A", "Alzheimer's Disease"),
    ("rs63750847", "A;A", "Alzheimer's Disease"),
    ("rs429358", "C;T", "Alzheimer's Disease"),
    ("rs145999145", "A;G", "Alzheimer's Disease"),
    ("rs63750847", "A;G", "Alzheimer's Disease"),
    
    ("rs1858830", "C;C", "Autism"),
    ("rs2710102", "C;C", "Autism"),
    ("rs7794745", "A;T", "Autism"),
    ("rs1322784", "C;C", "Autism"),
    ("rs1322784", "C;T", "Autism"),
    ("rs1322784", "T;T", "Autism"),
    ("rs265981", "A;G", "Autism"),
    ("rs4532", "C;T", "Autism"),
    ("rs686", "A;G", "Autism"),
    ("rs1143674", "A;A", "Autism"),
    ("rs6807362", "C;C", "Autism"),
    ("rs757972971", "A;A", "Autism"),
    ("rs2217262", "A;A", "Autism"),
    ("rs6766410", "A;A", "Autism"),
    ("rs6766410", "A;C", "Autism"),
    ("rs6766410", "C;C", "Autism"),
    ("rs1445442", "A;A", "Autism"),
    ("rs1445442", "A;G", "Autism"),
    ("rs1445442", "G;G", "Autism"),
    ("rs2421826", "C;T", "Autism"),
    ("rs2421826", "T;T", "Autism"),
    ("rs2421826", "C;C", "Autism"),
    ("rs1358054", "G;G", "Autism"),
    ("rs1358054", "G;T", "Autism"),
    ("rs1358054", "T;T", "Autism"),
    ("rs536861", "A;A", "Autism"),
    ("rs536861", "A;C", "Autism"),
    ("rs536861", "C;C", "Autism"),
    ("rs722628", "A;A", "Autism"),
    ("rs722628", "A;G", "Autism"),
    ("rs722628", "G;G", "Autism"),
    ("rs1858830", "C;G", "Autism"),
    ("rs2710102", "G;G", "Autism"),
    ("rs7794745", "T;T", "Autism"),
    ("rs265981", "G;G", "Autism"),
    ("rs4532", "T;T", "Autism"),
    ("rs686", "A;A", "Autism"),
    ("rs1143674", "A;A", "Autism"),
    ("rs757972971", "A;A", "Autism"),
    
    ("rs27388", "A;A", "Schizophrenia"),
    ("rs2270641", "G;G", "Schizophrenia"),
    ("rs4129148", "C;C", "Schizophrenia"),
    ("rs28694718", "A;A", "Schizophrenia"),
    ("rs6422441", "C;C", "Schizophrenia"),
    ("rs28414810", "C;C", "Schizophrenia"),
    ("rs6603272", "G;G", "Schizophrenia"),
    ("rs17883192", "C;C", "Schizophrenia"),
    ("rs165599", "G;G", "Schizophrenia"),
    ("rs27388", "A;G", "Schizophrenia"),
    ("rs4129148", "C;G", "Schizophrenia"),
    ("rs28694718", "A;G", "Schizophrenia"),
    ("rs6422441", "C;T", "Schizophrenia"),
    ("rs28414810", "C;G", "Schizophrenia"),
    ("rs6603272", "G;T", "Schizophrenia"),
    ("rs17883192", "C;G", "Schizophrenia"),
    
    ("rs3758391", "C;T", "Longevity"),
    ("rs5882", "A;A", "Longevity"),
    ("rs1042522", "C;C", "Longevity"),
    ("rs3803304", "C;C", "Longevity"),
    ("rs3803304", "C;G", "Longevity"),
    ("rs3803304", "G;G", "Longevity"),
    ("rs6873545", "C;C", "Longevity"),
    ("rs4590183", "C;C", "Longevity"),
    ("rs1556516", "C;C", "Longevity"),
    ("rs1556516", "C;G", "Longevity"),
    ("rs1556516", "G;G", "Longevity"),
    ("rs7137828", "C;C", "Longevity"),
    ("rs7137828", "C;T", "Longevity"),
    ("rs7137828", "T;T", "Longevity"),
    ("rs1627804", "C;C", "Longevity"),
    ("rs1627804", "A;A", "Longevity"),
    ("rs1627804", "A;C", "Longevity"),
    ("rs7844965", "A;G", "Longevity"),
    ("rs7844965", "A;A", "Longevity"),
    ("rs7844965", "G;G", "Longevity"),
    ("rs61978928", "C;C", "Longevity"),
    ("rs61978928", "C;T", "Longevity"),
    ("rs61978928", "T;T", "Longevity"),
    ("rs28926173", "C;C", "Longevity"),
    ("rs28926173", "C;T", "Longevity"),
    ("rs28926173", "T;T", "Longevity"),
    ("rs146254978", "C;C", "Longevity"),
    ("rs146254978", "C;T", "Longevity"),
    ("rs146254978", "T;T", "Longevity"),
    ("rs139137459", "A;G", "Longevity"),
    ("rs139137459", "A;A", "Longevity"),
    ("rs139137459", "G;G", "Longevity"),
    ("rs3758391", "T;T", "Longevity"),
    ("rs5882", "A;G", "Longevity"),
    ("rs1042522", "C;G", "Longevity"),
    
    ("rs333", ("46373456",), "Immunity"),
    
    ("rs28379706", "T;T", "Intelligence"),
    ("rs28379706", "C;T", "Intelligence"),
    ("rs28379706", "C;C", "Intelligence"),
    ("rs363039", "A;G", "Intelligence"),
    ("rs4680", "A;A", "Intelligence"),
    ("rs363039", "C;C", "Intelligence"),
    
    ("rs1815739", "C;C", "Muscular Performance"),
    ("rs1805086", "C;C", "Muscular Performance"),
    ("rs1815739", "C;T", "Muscular Performance"),
    ("rs1805086", "C;T", "Muscular Performance"),
    ("rs1815739", "T;T", "Muscular Performance"),
    
    ("rs4570625", "G;G", "OCD"),
    ("rs4565946", "C;C", "OCD"),
    
    ("rs1801131", "A;C", "Metabolism"),
    ("rs1801131", "C;C", "Metabolism"),
    ("rs1801133", "C;T", "Metabolism"),
    ("rs1801133", "T;T", "Metabolism"),
    ("rs2282679", "A;C", "Metabolism"),
    ("rs2282679", "C;C", "Metabolism"),
    ("rs12785878", "G;T", "Metabolism"),
    ("rs12785878", "T;T", "Metabolism"),
    ("rs1799945", "F;G", "Metabolism"),
    ("rs4988235", "C;C", "Metabolism"),
    ("rs182549", "C;C", "Metabolism"),
    ("rs2187668", "A;A", "Metabolism"),
    ("rs2187668", "A;G", "Metabolism"),
    ("rs5030858", "T;T", "Metabolism"),
    ("rs72921001", "C;C", "Metabolism"),
    ("rs7903146", "T;T", "Metabolism"),
    ("rs7903146", "C;C", "Metabolism"),
    ("rs7903146", "C;T", "Metabolism"),
    ("rs662799", "A;G", "Metabolism"),
    ("rs662799", "G;G", "Metabolism"),
    ("rs13119723", "A;A", "Metabolism"),
    ("rs13119723", "A;G", "Metabolism"),
    ("rs13119723", "G;G", "Metabolism"),
    ("rs6822844", "G;G", "Metabolism"),
    ("rs3184504", "C;T", "Metabolism"),
    ("rs3184504", "T;T", "Metabolism"),
    
    ("rs12913832", "A;A", "Eyes"),
    ("rs12913832", "A;G", "Eyes"),
    ("rs12913832", "G;G", "Eyes"),
    ("rs28938473", "T;T", "Eyes"),
    ("rs61753033", "T;T", "Eyes"),
    ("rs61753034", "T;T", "Eyes"),
    ("rs4778241", "A;A", "Eyes"),
    ("rs4778241", "A;C", "Eyes"),
    ("rs4778241", "C;C", "Eyes"),
    ("rs7495174", "A;A", "Eyes"),
    ("rs1129038", "A;A", "Eyes"),
    ("rs1129038", "A;G", "Eyes"),
    ("rs1129038", "G;G", "Eyes"),
    ("rs916977", "A;A", "Eyes"),
    ("rs916977", "A;G", "Eyes"),
    ("rs916977", "G;G", "Eyes"),
    ("rs1667394", "A;A", "Eyes"),
    
    ("rs6152", "A;A", "Hair"),
    ("rs6152", "A;G", "Hair"),
    ("rs6152", "A;", "Hair"),
    ("rs6152", "G;G", "Hair"),
    ("rs1805009", "C;C", "Hair"),
    ("rs1805009", "C;G", "Hair"),
    ("rs1805007", "C;T", "Hair"),
    ("rs1805007", "T;T", "Hair"),
    ("rs1805008", "C;T", "Hair"),
    ("rs1805008", "T;T", "Hair"),
    ("rs1805006", "A;A", "Hair"),
    ("rs1805006", "A;C", "Hair"),
    ("rs11547464", "A;A", "Hair"),
    ("rs11547464", "A;G", "Hair"),
    ("rs35264875", "T;T", "Hair"),
    ("rs7349332", "T;T", "Hair"),
    ("rs11803731", "T;T", "Hair"),
    ("rs17646946", "A;A", "Hair"),
    ("rs1667394", "A;A", "Hair"),
]

RSIDS = [Rsid(x[0], x[1], x[2]) for x in _RSIDS]

def parse_vcf(path: str) -> dict[str, tuple[str, str]]:
    """parse a VCF file into a dictionary with rsid as keys and (ref, alt) as values"""
    rsid_dict = {}
    with open(path) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue  # Skip headers
            cols = line.strip().split("\t")
            if cols[2].startswith("rs"):  # Check if it's an rsid
                rsid = cols[2]
                ref = cols[3]
                alt = cols[4]
                rsid_dict[rsid] = (ref, alt)
    return rsid_dict

def scan_genes(rsid_dict: dict[str, tuple[str, str]]) -> list[Rsid]:
    """scan the given gene dictionary, returning a list of RSIDS that were detected"""
    detected_rsids = []
    for rsid in RSIDS:
        if rsid.rsid in rsid_dict:
            ref, alt = rsid_dict[rsid.rsid]
            print(f"RSID {rsid.rsid} found in VCF file. Checking alleles...")
            if rsid.detect(ref, alt):
                detected_rsids.append(rsid)
    return detected_rsids

if __name__ == "__main__":
    filename = (
        sys.argv[1]
        if len(sys.argv) > 1
        else input("Enter the location of the VCF file: ")
    )
    rsid_dict = parse_vcf(filename)
    rsids = scan_genes(rsid_dict)
    if not rsids:
        print("No matching RSIDs found in the VCF file.")
    for rsid in rsids:
        print(rsid.info + "\n")
