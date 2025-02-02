# MiniClustV1

MiniClustV1 is a graphical user interface (GUI) for clustering BCR trees using the `miniclust` algorithm. This GUI allows users to select folders, specify output files, and visualize the resulting clusters and BCR trees.

## Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/your-username/miniclustv1.git
    cd miniclustv1
    ```

2. Install the required dependencies:
    ```bash
    pip install -r requirements.txt
    ```

3. Ensure you have the following dependencies installed:
    - Python 3.x
    - tkinter
    - pickle
    - Bio
    - biopython
    - ete3
    - matplotlib
    - mplcursors
    - numpy
    - pandas
    - scikit_learn
    - scipy
    - screeninfo

## Usage

1. Run the `MiniClustV1` GUI:
    ```bash
    python gui_miniclustV1.py
    ```

2. Interact with the GUI:
    - **Folder**: Select the folder containing the input BCR trees.
    - **Output Dismat File**: Specify the output dissimilarity matrix file.
    - **Output Cluster File**: Specify the output cluster file.
    - **B Parameter**: Set the value for the B parameter.

3. Execute `MiniClust`:
    - Click the "Run MiniClust" button to perform the clustering.
    - The optimal number of clusters will be displayed.

4. Visualize the Trees and Clusters:
    - **Show Trees**: Click the "Show Trees" button to display all trees.
    - **Number of Clusters**: Specify the number of clusters and click the "Show Clusters" button to visualize the clusters with `ete3`.


## Authors

- [Reza Hassanzadeh](https://github.com/Reza-HZ)


