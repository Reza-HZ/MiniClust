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
    - functionsV1 (containing miniclust, show_all_trees, and plot_with_ete3)

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

## Functions

- `run_miniclust()`: Runs the `miniclust` function and saves the output trees and node weights.
- `run_show_all_trees()`: Visualizes all trees using the specified parameters.
- `run_plot_with_ete3()`: Plots the clusters using the `ete3` library.
- `browse_folder()`: Opens a dialog to select a folder.
- `browse_file(entry_var)`: Opens a dialog to select a file for the given entry variable.

## Authors

- [Reza Hassanzadeh](https://github.com/Reza-HZ)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgements

- Special thanks to the developers of the `functionsV1` library for providing essential functions used in this project.

