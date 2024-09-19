
### **Datatypes and Operators**

#### **Introduction**

In Python, **data types** classify and determine the operations that can be performed on data. Understanding data types is crucial because they dictate how the interpreter will process and store information. **Operators** are symbols that perform operations on variables and values.

#### **Basic Data Types**

1. **Integers (`int`)**: Whole numbers, positive or negative, without decimals.

   ```python
   num_genes = 20000  # Approximate number of genes in humans
   ```

2. **Floating-point numbers (`float`)**: Numbers with decimal points.

   ```python
   expression_level = 12.5  # Arbitrary units
   ```

3. **Strings (`str`)**: Sequences of characters enclosed in quotes.

   ```python
   dna_sequence = "ATGCGATACGCTTACG"
   ```

4. **Booleans (`bool`)**: Logical values indicating `True` or `False`.

   ```python
   is_coding_sequence = True
   ```

5. **NoneType (`None`)**: Represents the absence of a value.

   ```python
   protein_structure = None  # No structure assigned yet
   ```

#### **Operators**

1. **Arithmetic Operators**

   - **Addition (`+`)**

     ```python
     total_mass = mass_protein1 + mass_protein2
     ```

   - **Subtraction (`-`)**

     ```python
     mass_difference = mass_protein1 - mass_protein2
     ```

   - **Multiplication (`*`)**

     ```python
     total_concentration = concentration * volume
     ```

   - **Division (`/`)**

     ```python
     molarity = moles / volume
     ```

   - **Exponentiation (`**`)**

     ```python
     area = radius ** 2 * 3.1416
     ```

   - **Modulus (`%`)**: Remainder of division.

     ```python
     remainder = 17 % 3  # Returns 2
     ```

   - **Floor Division (`//`)**: Division that rounds down to the nearest whole number.

     ```python
     quotient = 17 // 3  # Returns 5
     ```

2. **Comparison Operators**

   - **Equal to (`==`)**
   - **Not equal to (`!=`)**
   - **Greater than (`>`)**
   - **Less than (`<`)**
   - **Greater than or equal to (`>=`)**
   - **Less than or equal to (`<=`)**

   ```python
   if expression_level > threshold:
       print("Gene is overexpressed.")
   ```

3. **Logical Operators**

   - **Logical AND (`and`)**
   - **Logical OR (`or`)**
   - **Logical NOT (`not`)**

   ```python
   if is_coding_sequence and has_promoter:
       print("Sequence is likely a gene.")
   ```

#### **Examples in Biology**

**Calculating GC Content**

```python
dna_sequence = "ATGCGATACGCTTACG"

g_count = dna_sequence.count('G')
c_count = dna_sequence.count('C')
gc_content = ((g_count + c_count) / len(dna_sequence)) * 100

print(f"GC Content: {gc_content:.2f}%")
```

**Using Operators in Conditions**

```python
expression_level = 12.5  # Arbitrary units

if expression_level >= 10 and is_coding_sequence:
    print("Gene is significantly expressed.")
else:
    print("Gene expression is low or sequence is non-coding.")
```

---

### **Variables**

#### **Introduction**

Variables act as storage containers for data values. In Python, you assign a value to a variable using the equals sign (`=`). Unlike some languages, you don't need to declare the type of a variable explicitly.

#### **Variable Naming Rules**

- **Must start** with a letter or an underscore (`_`).
- **Cannot start** with a number.
- **Can contain** letters, numbers, and underscores.
- **Case-sensitive**: `gene` and `Gene` are different variables.

#### **Assigning Variables**

```python
gene_name = "TP53"
sequence_length = 1500  # Base pairs
```

#### **Multiple Assignments**

```python
x, y, z = 10, 20, 30
```

#### **Reassigning Variables**

Variables can be reassigned to new values or even different data types.

```python
variable = 5
variable = "Now I'm a string"
```

#### **Variables in Calculations**

```python
mass = 180  # grams
molar_mass = 180.16  # g/mol (for glucose)
moles = mass / molar_mass

print(f"Number of moles: {moles:.2f} mol")
```

#### **Best Practices**

- **Use descriptive names** for readability.

  ```python
  num_cells = 5000
  ```

- **Avoid using Python reserved words** as variable names.

  ```python
  # Avoid using names like 'for', 'if', 'else', 'def', etc.
  ```

- **Consistent naming conventions**: Use `snake_case` for variable names.

---

### **Strings**

#### **Introduction**

Strings are sequences of characters used to store text data. They can be enclosed in single quotes (`'...'`), double quotes (`"..."`), or triple quotes (`'''...'''` or `"""..."""` for multi-line strings).

#### **String Operations**

1. **Concatenation**

   ```python
   gene = "BRCA"
   number = "1"
   full_name = gene + number  # "BRCA1"
   ```

2. **Repetition**

   ```python
   separator = "-" * 10  # "----------"
   ```

3. **Indexing and Slicing**

   - **Indexing**

     ```python
     first_base = dna_sequence[0]  # 'A'
     ```

   - **Slicing**

     ```python
     first_codon = dna_sequence[0:3]  # 'ATG'
     ```

4. **String Methods**

   - **`upper()` and `lower()`**

     ```python
     sequence = "atgc"
     sequence = sequence.upper()  # "ATGC"
     ```

   - **`count()`**

     ```python
     num_adenines = dna_sequence.count('A')
     ```

   - **`find()`**

     ```python
     position = dna_sequence.find('CGT')  # Returns index or -1 if not found
     ```

   - **`replace()`**

     ```python
     rna_sequence = dna_sequence.replace('T', 'U')
     ```

#### **Formatting Strings**

- **Using `format()`**

  ```python
  message = "The gene {} has an expression level of {}.".format(gene_name, expression_level)
  ```

- **Using f-Strings (Python 3.6+)**

  ```python
  message = f"The gene {gene_name} has an expression level of {expression_level}."
  ```

#### **Examples in Biology**

**Complementing a DNA Sequence**

```python
def complement_dna(sequence):
    base_complements = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    complement = ''.join([base_complements[base] for base in sequence])
    return complement

dna_sequence = "ATGCGATACGCTTACG"
complement_sequence = complement_dna(dna_sequence)
print(f"Complementary Sequence: {complement_sequence}")
```

**Transcribing DNA to RNA**

```python
dna_sequence = "ATGCGATACGCTTACG"
rna_sequence = dna_sequence.replace('T', 'U')
print(f"RNA Sequence: {rna_sequence}")
```

---

### **Lists and Tuples**

#### **Lists**

Lists are ordered, mutable collections of items.

##### **Creating Lists**

```python
gene_list = ["BRCA1", "TP53", "EGFR"]
```

##### **Accessing Elements**

```python
first_gene = gene_list[0]  # "BRCA1"
```

##### **Modifying Lists**

- **Appending**

  ```python
  gene_list.append("MYC")
  ```

- **Inserting**

  ```python
  gene_list.insert(1, "PTEN")
  ```

- **Removing**

  ```python
  gene_list.remove("EGFR")
  ```

##### **Slicing**

```python
sublist = gene_list[1:3]
```

##### **List Comprehensions**

```python
squared_numbers = [x**2 for x in range(5)]  # [0, 1, 4, 9, 16]
```

##### **Examples in Biology**

**Creating a List of Nucleotide Counts**

```python
dna_sequence = "ATGCGATACGCTTACG"
nucleotide_counts = [dna_sequence.count(nuc) for nuc in ['A', 'T', 'G', 'C']]
print(nucleotide_counts)  # [4, 4, 4, 4]
```

#### **Tuples**

Tuples are ordered, immutable collections.

##### **Creating Tuples**

```python
codon = ('AUG', 'Methionine')
```

##### **Accessing Elements**

```python
start_codon = codon[0]  # "AUG"
```

##### **Immutability**

```python
codon[0] = 'UGA'  # This will raise a TypeError
```

##### **Examples in Biology**

**Storing Chromosome Locations**

```python
location = ('chr17', 43044295, 43125482)  # (chromosome, start, end)
```

---

### **Dictionary in Python**

#### **Introduction**

Dictionaries are unordered collections of key-value pairs. They are mutable and indexed by keys, which can be any immutable type.

#### **Creating Dictionaries**

```python
gene_expression = {
    "BRCA1": 75.0,
    "TP53": 60.5,
    "EGFR": 90.2
}
```

#### **Accessing Values**

```python
brca1_expression = gene_expression["BRCA1"]
```

#### **Adding or Updating Entries**

```python
gene_expression["MYC"] = 110.7  # Add new gene
gene_expression["BRCA1"] = 80.0  # Update existing gene
```

#### **Removing Entries**

```python
del gene_expression["TP53"]
```

#### **Dictionary Methods**

- **`keys()`**

  ```python
  genes = gene_expression.keys()
  ```

- **`values()`**

  ```python
  expressions = gene_expression.values()
  ```

- **`items()`**

  ```python
  for gene, expression in gene_expression.items():
      print(f"{gene}: {expression}")
  ```

#### **Examples in Biology**

**Codon to Amino Acid Mapping**

```python
codon_table = {
    'UUU': 'Phenylalanine',
    'UUC': 'Phenylalanine',
    'UUA': 'Leucine',
    'UUG': 'Leucine',
    # ... more codons
}

codon = 'AUG'
amino_acid = codon_table.get(codon, 'Unknown')
print(f"The codon {codon} codes for {amino_acid}.")
```

**Counting Amino Acids in a Protein**

```python
from collections import Counter

protein_sequence = "MVLSPADKTNVKAAW"
amino_acid_counts = Counter(protein_sequence)
print(amino_acid_counts)
```

---

### **Conditional Statements**

#### **Introduction**

Conditional statements allow you to execute code blocks based on certain conditions, controlling the flow of your program.

#### **`if`, `elif`, and `else`**

```python
if condition1:
    # Execute this block
elif condition2:
    # Execute this block
else:
    # Execute this block
```

#### **Examples in Biology**

**Assessing Gene Expression**

```python
expression_level = 85  # Arbitrary units

if expression_level > 100:
    print("High expression")
elif 50 < expression_level <= 100:
    print("Moderate expression")
else:
    print("Low expression")
```

**Analyzing Mutation Impact**

```python
mutation_type = "frameshift"

if mutation_type == "missense":
    print("Missense mutation: amino acid substitution.")
elif mutation_type == "nonsense":
    print("Nonsense mutation: introduces a stop codon.")
elif mutation_type == "frameshift":
    print("Frameshift mutation: alters the reading frame.")
else:
    print("Unknown mutation type.")
```

---

### **Loops in Python**

#### **Introduction**

Loops are used to execute a block of code repeatedly.

#### **`for` Loops**

Used for iterating over a sequence.

```python
for item in sequence:
    # Execute this block
```

##### **Example**

```python
genes = ["BRCA1", "TP53", "EGFR"]

for gene in genes:
    print(f"Analyzing {gene}")
```

#### **`while` Loops**

Executes as long as a condition is `True`.

```python
while condition:
    # Execute this block
```

##### **Example**

```python
count = 0

while count < 5:
    print(f"Count is {count}")
    count += 1
```

#### **Loop Control Statements**

- **`break`**: Exits the loop.

  ```python
  for num in range(10):
      if num == 5:
          break
      print(num)
  ```

- **`continue`**: Skips to the next iteration.

  ```python
  for num in range(10):
      if num % 2 == 0:
          continue
      print(num)  # Prints odd numbers
  ```

#### **Examples in Biology**

**Calculating Multiple GC Contents**

```python
dna_sequences = ["ATGCGATACGCTTACG", "GGCATCGTACGATCGA", "TTATCGCGATCGATGC"]

for seq in dna_sequences:
    gc_content = ((seq.count('G') + seq.count('C')) / len(seq)) * 100
    print(f"Sequence: {seq}, GC Content: {gc_content:.2f}%")
```

**Simulating Enzyme Kinetics**

```python
substrate_concentration = 0.1  # Starting concentration
km = 0.5
vmax = 1.2

while substrate_concentration <= 1.0:
    rate = (vmax * substrate_concentration) / (km + substrate_concentration)
    print(f"[S]: {substrate_concentration}, Rate: {rate:.2f}")
    substrate_concentration += 0.1
```

---

### **Functions**

#### **Introduction**

Functions are reusable blocks of code designed to perform a single, related action.

#### **Defining Functions**

```python
def function_name(parameters):
    # Function body
    return value
```

#### **Example**

```python
def calculate_molarity(moles, volume_liters):
    molarity = moles / volume_liters
    return molarity
```

#### **Calling Functions**

```python
moles = 0.5
volume = 1.0
molarity = calculate_molarity(moles, volume)
print(f"Molarity: {molarity} M")
```

#### **Default Parameters**

```python
def greet(name, message="Hello"):
    print(f"{message}, {name}!")

greet("Alice")  # Uses default message
greet("Bob", "Hi")  # Uses provided message
```

#### **Docstrings**

Provide documentation for your functions.

```python
def calculate_molarity(moles, volume_liters):
    """
    Calculates molarity given moles and volume in liters.

    Parameters:
    moles (float): Number of moles
    volume_liters (float): Volume in liters

    Returns:
    float: Molarity in moles per liter (M)
    """
    return moles / volume_liters
```

#### **Examples in Biology**

**Calculating Enzyme Activity**

```python
def enzyme_activity(vmax, substrate_concentration, km):
    """
    Calculates enzyme activity using the Michaelis-Menten equation.

    Parameters:
    vmax (float): Maximum rate
    substrate_concentration (float): Substrate concentration
    km (float): Michaelis constant

    Returns:
    float: Reaction rate
    """
    rate = (vmax * substrate_concentration) / (km + substrate_concentration)
    return rate

rate = enzyme_activity(vmax=1.2, substrate_concentration=0.5, km=0.3)
print(f"Enzyme activity: {rate:.2f}")
```

**Translating DNA to Protein**

```python
from Bio.Seq import Seq

def translate_sequence(dna_seq):
    """
    Translates a DNA sequence into a protein sequence.

    Parameters:
    dna_seq (str): DNA sequence

    Returns:
    str: Protein sequence
    """
    dna = Seq(dna_seq)
    protein = dna.translate(to_stop=True)
    return str(protein)

dna_sequence = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
protein_sequence = translate_sequence(dna_sequence)
print(f"Protein Sequence: {protein_sequence}")
```

---

### **Classes and Objects**

#### **Introduction**

Classes define the blueprint for objects, encapsulating data and functions that operate on data.

#### **Defining a Class**

```python
class Organism:
    def __init__(self, name, kingdom):
        self.name = name
        self.kingdom = kingdom

    def display_info(self):
        print(f"{self.name} belongs to the {self.kingdom} kingdom.")
```

#### **Creating an Object**

```python
human = Organism("Homo sapiens", "Animalia")
human.display_info()
```

#### **Inheritance**

```python
class Animal(Organism):
    def __init__(self, name, species):
        super().__init__(name, "Animalia")
        self.species = species

    def display_species(self):
        print(f"{self.name} is a {self.species}.")
```

```python
dog = Animal("Dog", "Canis lupus familiaris")
dog.display_info()
dog.display_species()
```

#### **Examples in Biology**

**Modeling Genes**

```python
class Gene:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def gc_content(self):
        gc = (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) * 100
        return gc

class ProteinCodingGene(Gene):
    def translate(self):
        from Bio.Seq import Seq
        dna_seq = Seq(self.sequence)
        protein_seq = dna_seq.translate(to_stop=True)
        return str(protein_seq)

brca1 = ProteinCodingGene("BRCA1", "ATGCGATACGCTTACG")
print(f"GC Content of {brca1.name}: {brca1.gc_content():.2f}%")
print(f"Protein Sequence: {brca1.translate()}")
```

---

### **File Handling in Python**

#### **Introduction**

File handling allows you to read from and write to files, essential for data processing.

#### **Reading a File**

```python
with open('data.txt', 'r') as file:
    content = file.read()
    print(content)
```

#### **Writing to a File**

```python
with open('results.txt', 'w') as file:
    file.write("Analysis Results\n")
    file.write("Gene Expression Levels\n")
```

#### **Appending to a File**

```python
with open('results.txt', 'a') as file:
    file.write("BRCA1: 75.0\n")
```

#### **Working with CSV Files**

```python
import csv

with open('gene_expression.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        print(row)
```

#### **Examples in Biology**

**Reading FASTA Files**

```python
from Bio import SeqIO

with open('sequences.fasta', 'r') as fasta_file:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        print(f"ID: {record.id}")
        print(f"Sequence: {record.seq}")
```

**Processing Experimental Data**

```python
import csv

with open('experiment_data.csv', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        sample_id = row['SampleID']
        value = float(row['Measurement'])
        print(f"Sample {sample_id}: {value}")
```

**Writing Analysis Results**

```python
results = {
    "Sample1": 0.95,
    "Sample2": 0.85,
    "Sample3": 0.75
}

with open('analysis_results.txt', 'w') as file:
    file.write("Sample\tValue\n")
    for sample, value in results.items():
        file.write(f"{sample}\t{value}\n")
```

---

### **Conclusion**

These foundational concepts in Python are essential tools for biologists engaging in computational work. By understanding data types, control structures, functions, and file handling, you can automate analyses, process large datasets, and model complex biological systems.

