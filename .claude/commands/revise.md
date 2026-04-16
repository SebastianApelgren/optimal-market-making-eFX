I want you to help me revise a chapter of my thesis based on revision notes.

Arguments: $ARGUMENTS

Parse the arguments as: <chapter-file> <notes-file> <group-letter>

For example: `report/chapters/ode.tex notes/ode_chapter_revisions.md C`

## Workflow

1. **Read both files.** Read the chapter file and the revision notes file in full.

2. **Extract the target group.** Find the group matching the given letter (e.g., "Group C") in the notes file. If the group is not found, list the available groups and ask which one I meant.

3. **Skip completed items.** Items marked `[x]` are already done. Only work on items marked `[ ]`.

4. **For each open item, present:**
   - **Location:** Where in the chapter this item lives. Use the section/subsection name (e.g., "Section 3.4 — Parameters") and quote the first few words of the affected sentence so it's easy to find. For example: `Section 3.4 (Parameters), sentence starting "we use the parameter set from..."`.
   - A short summary of the issue (1-2 sentences)
   - Your own proposed fix. The notes may contain suggested fixes, but treat those as a starting point only. Critically evaluate whether the suggested fix is correct and well-phrased, and propose your own version. Respect the writing style guidelines in CLAUDE.md.
   - Show the specific old text and new text so I can see exactly what would change.

5. **Wait for my approval.** Do NOT edit any files until I confirm which fixes to apply. I may want to modify your proposals or skip some items.

6. **Apply approved fixes.** After I approve, apply only the changes I accepted using the Edit tool.

## Important

- Respect the LaTeX style and the writing conventions from CLAUDE.md (no em dashes, no banned words, etc.).
- If a fix involves code changes (e.g., plotting scripts), identify the relevant source files and show proposed code changes alongside the LaTeX changes.
- If an item requires judgment about the economics or mathematics, flag your uncertainty rather than guessing.
