

#### **Running Jupyter Notebook**

**Introduction**

Jupyter Notebook is an open-source web application that allows you to create and share documents containing live code, equations, visualizations, and narrative text. It's an excellent tool for data analysis, visualization, and interactive exploration, making it highly suitable for biological research.

**Why Use Jupyter Notebook?**

- **Interactive Coding**: Write and execute code in small chunks (cells), making it easier to test and debug.
- **Data Visualization**: Integrate code with plots and charts inline.
- **Documentation**: Combine code with rich text, equations (using LaTeX), and images to document your analysis.
- **Reproducibility**: Share notebooks with colleagues to reproduce analyses.

**Launching Jupyter Notebook**

**Step 1: Activate Your Conda Environment**

If you have created a Conda environment (e.g., `bioinfo`), activate it:

```bash
conda activate bioinfo
```

**Step 2: Start Jupyter Notebook**

In the terminal or Anaconda Prompt, type:

```bash
jupyter notebook
```

- This command will start the Jupyter Notebook server and open the interface in your default web browser.
- The default URL is `http://localhost:8888/tree`.

**Navigating the Notebook Dashboard**

- **File Navigation**: The dashboard displays the contents of the current working directory.
- **Creating a New Notebook**: Click on **New** and select **Python 3** (or the version of Python you installed).
- **Managing Notebooks**: You can open, rename, move, or delete notebooks from the dashboard.

**Using Jupyter Notebook**

**1. Understanding the Interface**

- **Cells**: The basic unit of a notebook, which can contain code or Markdown text.
  - **Code Cells**: Execute Python code.
  - **Markdown Cells**: Contain formatted text, images, or equations.
- **Toolbar**: Contains buttons for common actions like saving, adding cells, and running code.
- **Menu Bar**: Offers options for file operations, editing, viewing, and more.

**2. Writing and Executing Code**

- **Code Execution**: Type your Python code in a cell and press **Shift + Enter** to execute.
- **Example**: Calculating the molecular weight of a DNA sequence.

  ```python
  from Bio.SeqUtils import molecular_weight
  from Bio.Seq import Seq

  dna_seq = Seq("ATGCGATACGCTTACG")
  mw = molecular_weight(dna_seq, seq_type='DNA')
  print(f"Molecular Weight: {mw} Da")
  ```

**3. Adding Markdown Cells**

- **Inserting Text**: Use Markdown cells to add explanations, headings, or images.
- **Formatting**: Utilize Markdown syntax for styling.

  - **Headers**:

    ```markdown
    # Header 1
    ## Header 2
    ### Header 3
    ```

  - **Bold and Italic**:

    ```markdown
    **Bold Text**
    *Italic Text*
    ```

  - **Lists**:

    ```markdown
    - Item 1
    - Item 2
    ```

  - **Equations**:

    ```markdown
    $$E = mc^2$$
    ```

**4. Visualization**

- **Plotting Data**: Use libraries like Matplotlib or Seaborn to create graphs.

  ```python
  import matplotlib.pyplot as plt

  gc_contents = [40, 50, 60, 55, 65]
  samples = ['Sample A', 'Sample B', 'Sample C', 'Sample D', 'Sample E']

  plt.bar(samples, gc_contents)
  plt.xlabel('Samples')
  plt.ylabel('GC Content (%)')
  plt.title('GC Content Across Samples')
  plt.show()
  ```

**5. Saving and Exporting Notebooks**

- **Auto-Save**: Jupyter automatically saves your notebook periodically.
- **Manual Save**: Click the save icon or press **Ctrl + S**.
- **Exporting**: Go to **File > Download As** to export the notebook in various formats (e.g., HTML, PDF).

**6. Managing Kernels**

- **What is a Kernel?**: The computational engine that executes the code contained in the notebook.
- **Restarting a Kernel**: If your notebook becomes unresponsive, you can restart the kernel via **Kernel > Restart**.
- **Shutting Down**: Close the notebook tab and shut down the kernel from the dashboard or via **File > Close and Halt**.

**Example: Analyzing a DNA Sequence**

Let's walk through a simple example of analyzing a DNA sequence for open reading frames (ORFs).

**Step 1: Import Libraries**

```python
from Bio.Seq import Seq
from Bio.SeqUtils import CodonUsage
```

**Step 2: Define the DNA Sequence**

```python
dna_seq = Seq("ATGCGATACGCTTACGTAGCTAGCTAGCTAA")
```

**Step 3: Find ORFs**

```python
def find_orfs(sequence):
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    for frame in range(3):
        trans = str(sequence[frame:].translate(to_stop=False))
        proteins = trans.split('*')
        for protein in proteins:
            if 'M' in protein:
                orfs.append(protein[protein.find('M'):])
    return orfs

orfs = find_orfs(dna_seq)
print("Open Reading Frames:")
for i, orf in enumerate(orfs):
    print(f"ORF {i+1}: {orf}")
```

**Step 4: Visualizing Codon Usage**

```python
from collections import Counter

codon_count = Counter([str(dna_seq[i:i+3]) for i in range(0, len(dna_seq)-2, 3)])
print("Codon Usage:")
for codon, count in codon_count.items():
    print(f"{codon}: {count}")
```

**Best Practices for Jupyter Notebooks**

- **Organize Your Notebook**: Use headings and sections to structure your content.
- **Document Your Code**: Add comments and Markdown cells to explain your analysis.
- **Keep Cells Small**: Break down your code into logical chunks for easier debugging.
- **Version Control**: Consider using Git for version control, though be cautious with large output data.

**Extensions and Customizations**

- **Installing Extensions**: Enhance functionality with Jupyter Notebook extensions.

  ```bash
  conda install -c conda-forge jupyter_contrib_nbextensions
  jupyter contrib nbextension install --user
  ```

- **Enabling Extensions**: Go to **Nbextensions** tab in the dashboard to enable desired extensions (e.g., Table of Contents, Codefolding).

**Collaborating with Others**

- **Sharing Notebooks**: You can share your `.ipynb` file directly.
- **Nbviewer**: Use [nbviewer](https://nbviewer.jupyter.org/) to share a static version of your notebook online.
- **GitHub Integration**: Notebooks can be rendered directly in GitHub repositories.

**Troubleshooting Common Issues**

- **Kernel Issues**: If code isn't executing, try restarting the kernel.
- **Missing Libraries**: Ensure you're working in the correct environment and that all required packages are installed.
- **Browser Compatibility**: Jupyter works best in Chrome or Firefox. If you encounter issues, try switching browsers.

**Additional Tools**

- **JupyterLab**: The next-generation user interface for Project Jupyter, offering all the familiar building blocks of the classic notebook in a more flexible and powerful user interface.

  - Launch JupyterLab:

    ```bash
    jupyter lab
    ```

- **Jupyter Widgets**: Add interactive widgets to your notebooks for enhanced interactivity.

  ```python
  import ipywidgets as widgets
  from IPython.display import display

  def f(x):
      return x

  widgets.interact(f, x=10);
  ```

**Conclusion**

Jupyter Notebook is a versatile tool that integrates code execution, data visualization, and documentation into a single, shareable document. It's particularly useful in biology for exploratory data analysis, sharing reproducible research, and collaborating with others.

