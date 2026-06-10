# bioinformatics-learning
My journey in bioinformatics, R, Linux and data analysis
---

## DAY 1 — Detailed Notes

---

### What is Terminal?
Terminal is the command center of my computer.
Instead of clicking, I type commands.
Every bioinformatics tool runs in Terminal.
To open: Command + Space → type Terminal → Enter

---

### What is VS Code?
VS Code is my code editor — where I write all my scripts.
It has the file explorer on the left, editor in the middle, terminal at the bottom.
I will write R scripts, Python scripts and bash scripts here.

---

### What is GitHub?
GitHub is where my code lives online.
My code exists in TWO places:
- Local = on my Mac (Desktop/bioinfo-learning)
- Remote = on GitHub (github.com/senaesen13/bioinformatics-learning)
Every day I push my work so it is saved online forever.

---

### What is a bash script?
A bash script is a file containing terminal commands.
Instead of typing commands one by one, I save them in a file and run them all at once.
In the future my RNA-seq pipeline will be a bash script.
File extension = .sh

---

### Terminal Commands I learned today

| Command | Meaning | Example |
|---|---|---|
| pwd | where am I right now? | pwd |
| ls | what is inside this folder? | ls |
| mkdir | create a new folder | mkdir week1 |
| cd | go into a folder | cd week1 |
| cd .. | go back one folder | cd .. |
| touch | create an empty file | touch practice.sh |
| bash | run a script | bash practice.sh |
| echo | print a message on screen | echo "Hello!" |

---

### Git Commands I learned today

| Command | Meaning |
|---|---|
| git clone | download a GitHub repo to my Mac |
| git add . | select all changed files to save |
| git commit -m "" | take a snapshot with a label |
| git push | send snapshot to GitHub online |

---

### Vocabulary

| Word | Meaning |
|---|---|
| repo | a folder that remembers every change, like a time machine |
| local | the copy on my Mac |
| remote | the copy on GitHub online |
| commit | a saved snapshot of my work |
| push | sending my work from Mac to GitHub |
| clone | downloading a GitHub repo to my Mac |
| token | a special password GitHub generates for terminal access |
| .sh | file extension for bash scripts |
| .md | file extension for markdown files |
---

## DAY 1 — Part 2 — New Terminal Commands

---

### cat — print entire file on screen
Command: cat genes.txt
What it does: reads the file and prints every line on screen
Real bioinformatics use: reading a FASTQ file or gene list
Example output:
BRCA1
TP53
EGFR
MYC

---

### head — show only the first lines
Command: head -3 genes.txt
What it does: shows only the first 3 lines
The number after - decides how many lines to show
Real bioinformatics use: some files have millions of lines
you cannot cat the whole thing, so you use head to peek at the first few lines
Example output:
BRCA1
TP53
EGFR

---

### tail — show only the last lines
Command: tail -3 genes.txt
What it does: shows only the last 3 lines
Real bioinformatics use: checking the end of a results file
Example output:
VHL
RB1
BRCA2

---

### grep — search for a word inside a file
Command: grep "BRCA" genes.txt
What it does: finds every line that contains the word BRCA
Real bioinformatics use: you have 50000 genes in a results file
you cannot scroll through all of them
grep finds your gene of interest instantly
Example output:
BRCA1
BRCA2

---

### pipe | — connect two commands together
Command: cat genes.txt | grep "BR"
What it does: reads the file AND searches at the same time
The | symbol sends the output of the left command into the right command
Think of it like a water pipe — data flows from left to right

Example of chaining 3 commands:
cat results.txt | grep "significant" | head -20
→ read results file
→ find only significant genes
→ show only top 20

---

### Why all of this matters
These commands work on any text file — gene lists, results, sequences
They are the building blocks of every bioinformatics pipeline
---

## DAY 1 — Part 3 — More Terminal Commands

---

### wc — count lines or words
Command: wc -l genes.txt
What it does: counts how many lines are in a file
-l = count lines
-w = count words
Example output: 10 genes.txt
Real bioinformatics use: FASTQ files have 4 lines per read
if wc -l shows 4000000 lines → you have 1000000 reads

---

### sort — sort lines alphabetically
Command: sort genes.txt
What it does: sorts all lines from A to Z
Example output: BRCA1, BRCA2, CDH1, EGFR, KRAS, MYC, PTEN, RB1, TP53, VHL
Real bioinformatics use: sorting gene lists, sorting results by name

---

### cp — copy a file
Command: cp genes.txt genes_backup.txt
What it does: creates an exact copy of the file
genes.txt = original file
genes_backup.txt = new copy
Real bioinformatics use: always make a backup before modifying important files

---

### mv — move or rename a file
Command: mv genes_backup.txt genes_copy.txt
What it does: renames the file
can also move a file to a different folder
Real bioinformatics use: organising and renaming output files

---

### rm — delete a file permanently
Command: rm genes_copy.txt
What it does: deletes the file forever
WARNING: no trash bin, no undo — be very careful
Real bioinformatics use: cleaning up large temporary files

---

### Combining commands with pipe
Example 1: sort then count
sort genes.txt | wc -l
→ sort the file, then count the lines

Example 2: find then count
grep "BRCA" genes.txt | wc -l
→ find all BRCA genes, then count how many

Example 3: three commands chained
cat results.txt | grep "significant" | wc -l
→ read file, find significant genes, count them

---

### Full command reference so far

| Command | Meaning |
|---|---|
| pwd | where am I right now? |
| ls | what is inside this folder? |
| mkdir | create a new folder |
| cd | go into a folder |
| cd .. | go back one folder |
| touch | create an empty file |
| cat | print entire file on screen |
| head -3 | show first 3 lines |
| tail -3 | show last 3 lines |
| grep | search for a word in a file |
| wc -l | count lines in a file |
| sort | sort lines alphabetically |
| cp | copy a file |
| mv | move or rename a file |
| rm | delete a file permanently |
| pipe = | connect two commands together |
| bash | run a script |