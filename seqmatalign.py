import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker
import matplotlib.backends.backend_pdf
import sys
import os
import traceback
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
            self.score_button = tk.Button(self.files_frame, text="Load score",
                                          command=lambda: self.loadFile("score", 3))
            self.score_button.grid(column=2, row=0)
            self.score_file = ttk.Label(self.files_frame)
            self.score_file.grid(column=2, row=1, padx=5)
            self.files_frame.grid(column=0, row=0)


        def loadFile(self, type, id):
            filetypes = None
            if(type=="sequence"):
                filetypes = "*.fna *.faa *.fasta"
            else:
                filetypes = "*.txt"

            try:
                file = filedialog.askopenfilename(title="Select "+type+" file",
                                                  filetypes=((type+" file", filetypes), ("all files","*.*")))
                if(len(file)):
                    with open(file, "r") as file_reader:
                        content = file_reader.readlines()
                        seq_name = ""
                        if(id==1 or id==2):
                            if(content[0].startswith(">")):
                                seq_name = content[0].strip()
                                content = "".join(content[1:]).strip().replace("\n", "")
                            else:
                                content = "".join(content).strip().replace("\n", "")

                            if(id==1):
                                self.seq1 = content
                                self.seq1_file.config(text=os.path.basename(file))
                                self.seq1_name.config(text=seq_name)
                            elif(id==2):
                                self.seq2 = content
                                self.seq2_file.config(text=os.path.basename(file))
                                self.seq2_name.config(text=seq_name)

                        elif(id==3):
                            self.scores = {}
                            letters = content[0].strip().split(" ")
                            first_index = 0
                            for line in content[1:]:
                                score = line.strip().split(" ")
                                self.scores[letters[first_index]] = {}
                                for second_index in range(len(score)):
                                    self.scores[letters[first_index]][letters[second_index]] = float(score[second_index])

                                first_index += 1

                            self.score_file.config(text=os.path.basename(file))



            except Exception:
                print(traceback.format_exc())
                messagebox.showerror("File error", "An error occurred while opening a file. The file is most likely just in a wrong/unknown format."
                                     " Check the file and try again or open a new file.")

        def backtrace(self, i, j, aligned_seq1, aligned_seq2):
            if((i == 0 and j == 0) or (self.args.algorithm == "sw" and self.matrix[i][j] == 0)):
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
                self.backtrace(i, j-1, ["-"]+aligned_seq1, [self.seq2[j]]+aligned_seq2)

            if("U" in self.backtrace_matrix[i][j]):
                try:
                    self.arrow_matrix[i][j]["U"].remove()
                except ValueError as ve:
                    pass

                plt.annotate(text="", xy=(j, i-0.15), xytext=(j, i-0.85), arrowprops=dict(arrowstyle="<-", color="b"))
                self.backtrace(i-1, j, [self.seq1[i]]+aligned_seq1, ["-"]+aligned_seq2)


        def highlight_cell(self, x, y, ax=None, **kwargs):
            rect = plt.Rectangle((x-.5, y-.5), 1,1, fill=False, **kwargs)
            ax = ax or plt.gca()
            ax.add_patch(rect)
            return rect


        def main(self):
            parser = argparse.ArgumentParser(description="Align two sequences with Needleman-Wunsch or Smith-Waterman and build the backtrace self.matrix")
            parser.add_argument("-s1", "--sequence1", help="Set the path to first sequence file (with or without header)", default="seq1.fasta")
            parser.add_argument("-s2", "--sequence2", help="Set the path to second sequence file (with or without header)", default="seq2.fasta")
            parser.add_argument("-alg", "--algorithm", help="Specify whether the Needleman-Wunsch (nw) or the Smith-Waterman (sw) algorithm should be used", choices=["nw", "sw"], default="nw")
            parser.add_argument("-f", "--function", help="Specify whether the score function should be minimized (distance, dis) or maximized (similarity, sim). If the Smith-Waterman algorithm is chosen, the function is always maximized", choices=["dis", "sim"], default="sim")
            parser.add_argument("-s", "--scores", help="Set the path to the score file", default="scores.txt")
            parser.add_argument("-g", "--gap", help="Specify the gap costs", default=-1, type=float)
            parser.add_argument("-ie", "--ignore-ends", help="Ignore overhanging ends in the alignment", action="store_true")
            self.args = parser.parse_args()

            self.seq1 = None
            with open(self.args.sequence1, "r") as seq_reader:
                content = seq_reader.readlines()
                if(content[0].startswith(">")):
                    self.seq1 = "".join(content[1:]).strip().replace("\n", "")
                else:
                    self.seq1 = "".join(content).strip().replace("\n", "")

            self.seq2 = None
            with open(self.args.sequence2, "r") as seq_reader:
                content = seq_reader.readlines()
                if(content[0].startswith(">")):
                    self.seq2 = "".join(content[1:]).strip().replace("\n", "")
                else:
                    self.seq2 = "".join(content).strip().replace("\n", "")

            self.scores = {}
            with open(self.args.scores, "r") as score_reader:
                content = score_reader.readlines()
                letters = content[0].strip().split(" ")
                first_index = 0
                for line in content[1:]:
                    score = line.strip().split(" ")
                    scores[letters[first_index]] = {}
                    for second_index in range(len(score)):
                        scores[letters[first_index]][letters[second_index]] = float(score[second_index])

                    first_index += 1

            alg_type = "global"
            if(self.args.algorithm == "sw"):
                self.args.function = "sim"
                alg_type = "local"

            self.seq1 = "-" + self.seq1
            self.seq2 = "-" + self.seq2
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
                if(self.args.algorithm == "nw"):
                    if(self.args.ignore_ends):
                        self.matrix[0][j] = 0.0
                    else:
                        self.matrix[0][j] = float(self.args.gap) * j

                    self.backtrace_matrix[0][j] = "L"
                    arrow = plt.annotate(text="", xy=(j-0.2, 0), xytext=(j-0.8, 0), arrowprops=dict(arrowstyle="<-", color="r"))
                    self.arrow_matrix[0][j]["L"] = arrow

                if(self.matrix[0][j].is_integer()):
                    ax.text(j, 0, str(int(self.matrix[0][j])), va='center', ha='center')
                else:
                    ax.text(j, 0, str(self.matrix[0][j]), va='center', ha='center')

                self.highlight_cell(j, 0, color="black", linewidth=1)

            for i in range(1, len(self.seq1)):
                if(self.args.algorithm == "nw"):
                    if(self.args.ignore_ends):
                        self.matrix[i][0] = 0.0
                    else:
                        self.matrix[i][0] = float(self.args.gap) * i

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

                    diagonal_value = self.matrix[i-1][j-1] + scores[self.seq1[i]][self.seq2[j]]
                    left_value = self.matrix[i][j-1] + self.args.gap
                    if(self.args.ignore_ends and i==len(self.seq1)-1):
                        left_value = self.matrix[i][j-1]

                    up_value = self.matrix[i-1][j] + self.args.gap
                    if(self.args.ignore_ends and j==len(self.seq2)-1):
                        up_value = self.matrix[i-1][j]

                    if(self.args.function == "sim"):
                        self.matrix[i][j] = max(diagonal_value, max(left_value, up_value))
                        if(self.args.algorithm == "sw"):
                            self.matrix[i][j] = max(0.0, self.matrix[i][j])
                    elif(self.args.function == "dis"):
                        self.matrix[i][j] = min(diagonal_value, min(left_value, up_value))

                    if(self.args.algorithm == "nw" or (self.args.algorithm == "sw" and self.matrix[i][j] > 0)):
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
            if(self.args.algorithm == "nw"):
                self.highlight_cell(len(self.seq2)-1, len(self.seq1)-1, color="green", linewidth=2)
                self.backtrace(len(self.seq1)-1, len(self.seq2)-1, [], [])
            elif((self.args.algorithm == "sw" and maximum_value != 0)):
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
            plt.savefig(alg_type+"_alignments.pdf", bbox_inches="tight")

    try:
        root = tk.Tk()
        root.style = ttk.Style()
        root.style.theme_use("clam")
        Application(root)
        root.mainloop()
    except Exception:
        print(traceback.format_exc())
