# ActiveGrid

GitHub repository created by Nathan Wei and Kevin Griffin, Summer 2015

Welcome to the ActiveGrid repository on GitHub!

Compiling and Running Code (from PaddleCode folder):
1) make cleaner - this command removes temporary and localized files, so they don't cause merge conflicts on the server
2) source source.txt - this runs make files for both the PaddleCode and lib folders
3) ./menuII - this runs the program. Press ctrl-C to exit at any time.

master is the main branch for this project. Substantial changes to the project should be made in this manner:
1) create a new branch off master
2) make changes in a series of commits to this branch
3) checkout and pull down from master (if anyone else has modified it while you were on your local branch)
4) go to your branch and merge master into it (this should take care of most merge conflicts; you may have to remove some residual ones yourself afterwards)
5) go to the repository website on GitHub and create a pull request for your branch
6) merge the pull request (and delete the branch if desired)

Here are some useful GitHub commands:
- git init (initialize repository)
- git checkout -b branchname (make new branch off master)
- git add filename (add file to commit)
- git commit -m "message" (commit changes you've made to the files you've added)
- git push origin branchname (push local commits/changes to the server)
- git checkout master (switch to master branch)
- git pull (pull down changes from master)
- git checkout branchname, git merge master (merge changes from master into your branch)
- git status (see information on the state of your local repository)

If you have any questions about architecture, implementations, or this method of version control, feel free to email us at nwei@princeton.edu or kevinpg@princeton.edu!