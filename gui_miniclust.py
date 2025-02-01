import tkinter as tk
from tkinter import filedialog, messagebox
import pickle
from functions import miniclust

def run_miniclust():
    folder = folder_var.get()
    folder += '/'
    output_dismat_file = dismat_file_var.get()
    output_cluster_file = cluster_file_var.get()
    num_clusters = int(num_clusters_var.get())
    show_clusters = show_clusters_var.get() == '1'
    show_simp_trees = show_simp_trees_var.get() == '1'
    b = int(b_var.get())

    if not folder or not output_dismat_file or not output_cluster_file:
        messagebox.showerror("Input Error", "Please fill in all required fields.")
        return

    try:
        new_trees, new_trees_nodes_weights = miniclust(
            folder,
            output_dismat_file,
            output_cluster_file,
            num_clusters=num_clusters,
            show_clusters=show_clusters,
            show_simp_trees=show_simp_trees,
            b=b
        )

        with open('output.pkl', 'wb') as file:
            pickle.dump((new_trees, new_trees_nodes_weights), file)
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

def browse_folder():
    folder = filedialog.askdirectory()
    if folder:
        folder_var.set(folder)

def browse_file(entry_var):
    file = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text Files", "*.txt")])
    if file:
        entry_var.set(file)

# Create the main window
root = tk.Tk()
root.title("MiniClustV0")
root.geometry("450x250")

# Variables
folder_var = tk.StringVar()
dismat_file_var = tk.StringVar()
cluster_file_var = tk.StringVar()
num_clusters_var = tk.StringVar(value="5")
show_clusters_var = tk.StringVar(value="0")
show_simp_trees_var = tk.StringVar(value="0")
b_var = tk.StringVar(value="0")

# Layout
frame = tk.Frame(root, padx=10, pady=10)
frame.pack(fill="both", expand=True)

# Folder selection
folder_label = tk.Label(frame, text="Folder:")
folder_label.grid(row=0, column=0, sticky="w")
folder_entry = tk.Entry(frame, textvariable=folder_var, width=40)
folder_entry.grid(row=0, column=1)
browse_folder_button = tk.Button(frame, text="Browse", command=browse_folder)
browse_folder_button.grid(row=0, column=2)

# Output dismat file
dismat_label = tk.Label(frame, text="Output Dismat File:")
dismat_label.grid(row=1, column=0, sticky="w")
dismat_entry = tk.Entry(frame, textvariable=dismat_file_var, width=40)
dismat_entry.grid(row=1, column=1)
dismat_button = tk.Button(frame, text="Browse", command=lambda: browse_file(dismat_file_var))
dismat_button.grid(row=1, column=2)

# Output cluster file
cluster_label = tk.Label(frame, text="Output Cluster File:")
cluster_label.grid(row=2, column=0, sticky="w")
cluster_entry = tk.Entry(frame, textvariable=cluster_file_var, width=40)
cluster_entry.grid(row=2, column=1)
cluster_button = tk.Button(frame, text="Browse", command=lambda: browse_file(cluster_file_var))
cluster_button.grid(row=2, column=2)

# Number of clusters
num_clusters_label = tk.Label(frame, text="Number of Clusters:")
num_clusters_label.grid(row=3, column=0, sticky="w")
num_clusters_entry = tk.Entry(frame, textvariable=num_clusters_var, width=10)
num_clusters_entry.grid(row=3, column=1, sticky="w")

# Show clusters
show_clusters_label = tk.Label(frame, text="Show Clusters:")
show_clusters_label.grid(row=4, column=0, sticky="w")
show_clusters_check = tk.Checkbutton(frame, variable=show_clusters_var, onvalue="1", offvalue="0")
show_clusters_check.grid(row=4, column=1, sticky="w")

# Show simple trees
show_simp_trees_label = tk.Label(frame, text="Show Minified Trees:")
show_simp_trees_label.grid(row=5, column=0, sticky="w")
show_simp_trees_check = tk.Checkbutton(frame, variable=show_simp_trees_var, onvalue="1", offvalue="0")
show_simp_trees_check.grid(row=5, column=1, sticky="w")

# B parameter
b_label = tk.Label(frame, text="B Parameter:")
b_label.grid(row=6, column=0, sticky="w")
b_entry = tk.Entry(frame, textvariable=b_var, width=10)
b_entry.grid(row=6, column=1, sticky="w")

# Run button
run_button = tk.Button(frame, text="Run", width=20, height=2, command=run_miniclust)
run_button.grid(row=7, column=0, columnspan=3, pady=10)

# Start the GUI event loop
root.mainloop()
