import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker
import matplotlib.backends.backend_pdf
import sys
import os
import traceback
import re
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

if __name__ == "__main__":

    class Application(tk.Frame):

        def __init__(self, master=None):
            tk.Frame.__init__(self, master)
            self.master = master
            self.initWindow()

        def initWindow(self):
            self.master.title("SeqMatAlign")
            self.pack(fill="both", expand=1)

            ####################################################################
            self.files_frame = tk.LabelFrame(self, text="Files", borderwidth=2, relief="groove")

            self.seq1 = None
            self.seq1_button = tk.Button(self.files_frame, text="Load first sequence",
                                         command=lambda: self.loadFile("sequence", 1))
            self.seq1_button.grid(column=0, row=0, padx=5)
            self.seq1_file = ttk.Label(self.files_frame)
            self.seq1_file.grid(column=0, row=1, padx=5, pady=5)
            self.seq1_name = ttk.Label(self.files_frame)
            self.seq1_name.grid(column=0, row=2, padx=5, pady=5)

            self.seq2 = None
            self.seq2_button = tk.Button(self.files_frame, text="Load second sequence",
                                         command=lambda: self.loadFile("sequence", 2))
            self.seq2_button.grid(column=1, row=0, padx=5)
            self.seq2_file = ttk.Label(self.files_frame)
            self.seq2_file.grid(column=1, row=1, padx=5, pady=5)
            self.seq2_name = ttk.Label(self.files_frame)
            self.seq2_name.grid(column=1, row=2, padx=5, pady=5)

            self.scores = None
            self.score_button = tk.Button(self.files_frame, text="Load score matrix",
                                          command=lambda: self.loadFile("score", 3))
            self.score_button.grid(column=2, row=0)
            self.score_file = ttk.Label(self.files_frame)
            self.score_file.grid(column=2, row=1)
            self.remove_score_button = tk.Button(self.files_frame, text="Remove score matrix",
                                                 command=lambda: self.removeScore())
            self.remove_score_button.grid(column=2, row=2)

            self.files_frame.grid(column=0, row=0)

            ####################################################################
            self.options_frame = tk.LabelFrame(self, text="Options", borderwidth=2, relief="groove")

            self.match_label = ttk.Label(self.options_frame, text="Match score:")
            self.match_label.grid(column=0, row=0)
            self.match_entry = tk.Entry(self.options_frame, width=4)
            self.match_entry.insert(0, 2)
            self.match_entry.grid(column=1, row=0, padx=5)

            self.mismatch_label = ttk.Label(self.options_frame, text="Mismatch score:")
            self.mismatch_label.grid(column=2, row=0)
            self.mismatch_entry = tk.Entry(self.options_frame, width=4)
            self.mismatch_entry.insert(0, 0)
            self.mismatch_entry.grid(column=3, row=0, padx=4)

            self.gap_label = ttk.Label(self.options_frame, text="Gap score:")
            self.gap_label.grid(column=4, row=0)
            self.gap_entry = tk.Entry(self.options_frame, width=4)
            self.gap_entry.insert(0, -1)
            self.gap_entry.grid(column=5, row=0)

            self.algo_label = ttk.Label(self.options_frame, text="Algorithm:")
            self.algo_label.grid(column=0, row=1, pady=5)
            self.algo_type = tk.StringVar()
            self.algo_option = ttk.OptionMenu(self.options_frame, self.algo_type, "NW", *["NW", "SW"])
            self.algo_option.grid(column=1, row=1, padx=4, pady=5)

            self.function_label = ttk.Label(self.options_frame, text="Function:")
            self.function_label.grid(column=2, row=1, pady=5)
            self.function_type = tk.StringVar()
            self.function_option = ttk.OptionMenu(self.options_frame, self.function_type, "Similarity", *["Similarity", "Distance"])
            self.function_option.grid(column=3, row=1, padx=4, pady=5)

            self.ignore_ends_var = tk.IntVar()
            self.ignore_ends = tk.Checkbutton(self.options_frame, text="Ignore ends", variable=self.ignore_ends_var)
            self.ignore_ends.grid(column=4, row=1, pady=5)

            self.options_frame.grid(column=0, row=1, pady=10)

            ####################################################################
            self.start_frame = tk.LabelFrame(self, text="SeqMatAlign", borderwidth=2, relief="groove")

            self.start_button = tk.Button(self.start_frame, text="Start aligning",
                                          command=lambda: self.align())
            self.start_button.grid(column=0, row=0, padx=10)

            self.start_frame.grid(column=0, row=2, pady=10)


        def loadFile(self, type, id):
            filetypes = None
            if(type=="sequence"):
                filetypes = "*.fna *.faa *.fasta"
            else:
                filetypes = "*.txt"

            file = filedialog.askopenfilename(title="Select "+type+" file",
                                              filetypes=((type+" file", filetypes), ("all files","*.*")))
            if(len(file)):
                try:
                    with open(file, "r") as file_reader:
                        content = file_reader.readlines()
                        seq = []
                        seq_name = "no header found"
                        if(id==1 or id==2):
                            for line in content:
                                line = line.strip()
                                if(line.startswith(">")):
                                    if(len(seq)):
                                        break

                                    seq_name = line

                                elif(not line.startswith("#")):
                                    seq.append(line)

                            if(id==1):
                                self.seq1 = "*" + "".join(seq).upper()
                                self.seq1_file.config(text=os.path.basename(file))
                                self.seq1_name.config(text=seq_name)
                            elif(id==2):
                                self.seq2 = "*" + "".join(seq).upper()
                                self.seq2_file.config(text=os.path.basename(file))
                                self.seq2_name.config(text=seq_name)

                        elif(id==3):
                            self.scores = {}
                            col_names = []
                            row_number = 0
                            first_row = True
                            for line in content:
                                line = line.strip()
                                if(not line.startswith("#")):
                                    if(first_row):
                                        col_names = list(filter(None, line.split(" ")))
                                        first_row = False
                                    else:
                                        row = list(filter(None, line.split(" ")))
                                        if(not row[0] == col_names[row_number]):
                                            messagebox.showerror("Invalid matrix format", "The row name doesn't match the column name at position {}.".format(row_number))
                                            return

                                        self.scores[row[0]] = {}
                                        for c in range(1, len(row)):
                                            if(not self.validateEntry(row[c])):
                                                messagebox.showerror("Invalid matrix entry", "Invalid number at row {} and column {}.".format(row_number, c-1))
                                                return

                                            self.scores[row[0]][col_names[c-1]] = float(row[c])

                                        row_number += 1

                            self.score_file.config(text=os.path.basename(file))
                            self.match_entry.config(state="disabled")
                            self.mismatch_entry.config(state="disabled")
                            self.gap_entry.config(state="disabled")
                except FileNotFoundError:
                    messagebox.showerror("File not accessible", "The file was either not found or unable to be accessed.")


        def removeScore(self):
            self.score_file["text"] = ""
            self.match_entry.config(state="normal")
            self.mismatch_entry.config(state="normal")
            self.gap_entry.config(state="normal")


        def validateEntry(self, entry):
            try:
                float(entry)
            except ValueError:
                print(traceback.format_exc())
                return False

            return True


        def backtrace(self, i, j, aligned_seq1, aligned_seq2):
            if((i == 0 and j == 0) or (self.algo_type.get() == "SW" and self.matrix[i][j] == 0)):
                self.alignments.append([aligned_seq1, aligned_seq2])
                return

            if("D" in self.backtrace_matrix[i][j]):
                try:
                    self.arrow_matrix[i][j]["D"].remove()
                except ValueError as ve:
                    pass

                plt.annotate(text="", xy=(j-0.2, i-0.15), xytext=(j-0.8, i-0.85), arrowprops=dict(arrowstyle="<-", color="b"))
                self.backtrace(i-1, j-1, [self.seq1[i]]+aligned_seq1, [self.seq2[j]]+aligned_seq2)

            if("L" in self.backtrace_matrix[i][j]):
                try:
                    self.arrow_matrix[i][j]["L"].remove()
                except ValueError as ve:
                    pass

                plt.annotate(text="", xy=(j-0.2, i), xytext=(j-0.8, i), arrowprops=dict(arrowstyle="<-", color="b"))
                self.backtrace(i, j-1, ["*"]+aligned_seq1, [self.seq2[j]]+aligned_seq2)

            if("U" in self.backtrace_matrix[i][j]):
                try:
                    self.arrow_matrix[i][j]["U"].remove()
                except ValueError as ve:
                    pass

                plt.annotate(text="", xy=(j, i-0.15), xytext=(j, i-0.85), arrowprops=dict(arrowstyle="<-", color="b"))
                self.backtrace(i-1, j, [self.seq1[i]]+aligned_seq1, ["*"]+aligned_seq2)


        def highlight_cell(self, x, y, ax=None, **kwargs):
            rect = plt.Rectangle((x-.5, y-.5), 1,1, fill=False, **kwargs)
            ax = ax or plt.gca()
            ax.add_patch(rect)
            return rect


        def align(self):
            try:
                if(not len(self.seq1_file["text"])):
                    messagebox.showerror("No first sequence file", "No file was given for the first sequence.")
                    return

                if(not len(self.seq2_file["text"])):
                    messagebox.showerror("No second sequence file", "No file was given for the second sequence.")
                    return

                if(not len(self.score_file["text"])):
                    if(not self.validateEntry(self.match_entry.get())):
                        messagebox.showerror("Invalid match cost", "Either the given match cost was not a number or the field was left empty.")
                        return

                    if(not self.validateEntry(self.mismatch_entry.get())):
                        messagebox.showerror("Invalid mismatch cost", "Either the given mismatch cost was not a number or the field was left empty.")
                        return

                    if(not self.validateEntry(self.gap_entry.get())):
                        messagebox.showerror("Invalid gap cost", "Either the given gap cost was not a number or the field was left empty.")
                        return

                    self.scores = {}
                    letters = set(self.seq1+self.seq2)
                    for letter1 in letters:
                        self.scores[letter1] = {}
                        for letter2 in letters:
                            if(letter1 == "*" or letter2 == "*"):
                                self.scores[letter1][letter2] = float(self.gap_entry.get())
                            elif(letter1==letter2):
                                self.scores[letter1][letter2] = float(self.match_entry.get())
                            else:
                                self.scores[letter1][letter2] = float(self.mismatch_entry.get())

                if(self.algo_type.get()=="SW"):
                    self.function_type.set("Similarity")

                self.matrix = [[0.0 for i in range(len(self.seq2))] for i in range(len(self.seq1))]
                self.backtrace_matrix = [["" for i in range(len(self.seq2))] for i in range(len(self.seq1))]
                self.arrow_matrix = [[{} for i in range(len(self.seq2))] for i in range(len(self.seq1))]
                maximum_value = -(sys.maxsize - 1)
                maximum_indices = []
                figsize = max(len(self.seq1), len(self.seq2))/2
                fig, ax = plt.subplots(figsize=(figsize, figsize))
                if(self.matrix[0][0].is_integer()):
                    ax.text(0, 0, str(int(self.matrix[0][0])), va='center', ha='center')
                else:
                    ax.text(0, 0, str(self.matrix[0][0]), va='center', ha='center')

                for j in range(1, len(self.seq2)):
                    if(self.algo_type.get() == "NW"):
                        if(self.ignore_ends_var.get()):
                            self.matrix[0][j] = 0.0
                        else:
                            self.matrix[0][j] = self.scores["*"][self.seq2[j]] * j

                        self.backtrace_matrix[0][j] = "L"
                        arrow = plt.annotate(text="", xy=(j-0.2, 0), xytext=(j-0.8, 0), arrowprops=dict(arrowstyle="<-", color="r"))
                        self.arrow_matrix[0][j]["L"] = arrow

                    if(self.matrix[0][j].is_integer()):
                        ax.text(j, 0, str(int(self.matrix[0][j])), va='center', ha='center')
                    else:
                        ax.text(j, 0, str(self.matrix[0][j]), va='center', ha='center')

                    self.highlight_cell(j, 0, color="black", linewidth=1)

                for i in range(1, len(self.seq1)):
                    if(self.algo_type.get() == "NW"):
                        if(self.ignore_ends_var.get()):
                            self.matrix[i][0] = 0.0
                        else:
                            self.matrix[i][0] = self.scores[self.seq1[i]]["*"] * i

                        self.backtrace_matrix[i][0] = "U"
                        arrow = plt.annotate(text="", xy=(0, i-0.15), xytext=(0, i-0.85), arrowprops=dict(arrowstyle="<-", color="r"))
                        self.arrow_matrix[i][0]["U"] = arrow

                    if(self.matrix[i][0].is_integer()):
                        ax.text(0, i, str(int(self.matrix[i][0])), va='center', ha='center')
                    else:
                        ax.text(0, i, str(self.matrix[i][0]), va='center', ha='center')

                    self.highlight_cell(0, i, color="black", linewidth=1)

                for i in range(1, len(self.seq1)):
                    for j in range(1, len(self.seq2)):
                        diagonal_value = self.matrix[i-1][j-1] + self.scores[self.seq1[i]][self.seq2[j]]
                        left_value = self.matrix[i][j-1] + self.scores["*"][self.seq2[j]]
                        if(self.ignore_ends_var.get() and i==len(self.seq1)-1):
                            left_value = self.matrix[i][j-1]

                        up_value = self.matrix[i-1][j] + self.scores[self.seq1[i]]["*"]
                        if(self.ignore_ends_var.get() and j==len(self.seq2)-1):
                            up_value = self.matrix[i-1][j]

                        if(self.function_type.get() == "Similarity"):
                            self.matrix[i][j] = max(diagonal_value, max(left_value, up_value))
                            if(self.algo_type.get() == "SW"):
                                self.matrix[i][j] = max(0.0, self.matrix[i][j])
                        elif(self.function_type.get() == "Distance"):
                            self.matrix[i][j] = min(diagonal_value, min(left_value, up_value))

                        if(self.algo_type.get() == "NW" or (self.algo_type.get() == "SW" and self.matrix[i][j] > 0)):
                            if(self.matrix[i][j] == diagonal_value):
                                self.backtrace_matrix[i][j] += "D"
                                arrow = plt.annotate(text="", xy=(j-0.2, i-0.15), xytext=(j-0.8, i-0.85), arrowprops=dict(arrowstyle="<-", color="r"))
                                self.arrow_matrix[i][j]["D"] = arrow

                            if(self.matrix[i][j] == left_value):
                                self.backtrace_matrix[i][j] += "L"
                                arrow = plt.annotate(text="", xy=(j-0.2, i), xytext=(j-0.8, i), arrowprops=dict(arrowstyle="<-", color="r"))
                                self.arrow_matrix[i][j]["L"] = arrow

                            if(self.matrix[i][j] == up_value):
                                self.backtrace_matrix[i][j] += "U"
                                arrow = plt.annotate(text="", xy=(j, i-0.15), xytext=(j, i-0.85), arrowprops=dict(arrowstyle="<-", color="r"))
                                self.arrow_matrix[i][j]["U"] = arrow

                            if(self.matrix[i][j] >= maximum_value):
                                if(self.matrix[i][j] > maximum_value):
                                    del maximum_indices[:]

                                maximum_value = self.matrix[i][j]
                                maximum_indices.append([i, j])

                        if(self.matrix[i][j].is_integer()):
                            ax.text(j, i, str(int(self.matrix[i][j])), va='center', ha='center')
                        else:
                            ax.text(j, i, str(self.matrix[i][j]), va='center', ha='center')

                        self.highlight_cell(j, i, color="black", linewidth=1)

                ax.matshow(self.matrix, aspect='auto', cmap=ListedColormap(["w"]))
                ax.set_yticklabels([""]+list(self.seq1))
                ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
                ax.set_xticklabels([""]+list(self.seq2))
                ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
                self.alignments = []
                if(self.algo_type.get() == "NW"):
                    self.highlight_cell(len(self.seq2)-1, len(self.seq1)-1, color="green", linewidth=2)
                    self.backtrace(len(self.seq1)-1, len(self.seq2)-1, [], [])
                elif((self.algo_type.get() == "SW" and maximum_value != 0)):
                    for indices in maximum_indices:
                        self.highlight_cell(indices[1], indices[0], color="green", linewidth=2)
                        self.backtrace(indices[0], indices[1], [], [])


                aligments_str = ""
                for alignment in self.alignments:
                    if(not len(aligments_str)):
                        aligments_str = str(alignment[0]) + "\n" + str(alignment[1])
                    else:
                        aligments_str += "\n\n" + str(alignment[0]) + "\n" + str(alignment[1])

                plt.gcf().text(0.93, 0.5, aligments_str, verticalalignment="center", bbox=dict(boxstyle='round', facecolor='white', alpha=0.15))
                plt.savefig("../../"+self.algo_type.get()+"_alignments.pdf", bbox_inches="tight")
            except Exception:
                messagebox.showerror("Unexpected error", "Something unexpected happened. Please ask the admin.")



    try:
        root = tk.Tk()
        def on_closing():
            root.quit()
            root.destroy()

        root.style = ttk.Style()
        root.style.theme_use("clam")
        Application(root)
        root.protocol("WM_DELETE_WINDOW", on_closing)
        root.mainloop()
    except Exception:
        messagebox.showerror("Unexpected error", "Something unexpected happened. Please ask the admin.")
