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