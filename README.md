# ActiveGrid

###### GitHub repository created by Nathan Wei and Kevin Griffin, Summer 2015

##Welcome to the ActiveGrid repository on GitHub!

#### Compiling and Running Code (from PaddleCode folder):
1. make cleaner - _this command removes temporary and localized files, so they don't cause merge conflicts on the server_
2. source source.txt - _this runs make files for both the PaddleCode and lib folders_
3. ./menuII - _this runs the program. Press ctrl-C to exit at any time._

##### master is the main branch for this project. Substantial changes to the project should be made in this manner:
1. create a new branch off master
2. make changes in a series of commits to this branch
3. checkout and pull down from master (if anyone else has modified it while you were on your local branch)
4. go to your branch and merge master into it (this should take care of most merge conflicts; you may have to remove some residual ones yourself afterwards)
5. go to the repository website on GitHub and create a pull request for your branch
6. merge the pull request (and delete the branch if desired)

##### Here are some useful GitHub commands:
- git init _(initialize repository)_
- git checkout -b branchname _(make new branch off master)_
- git add filename _(add file to commit)_
- git commit -m "message" _(commit changes you've made to the files you've added)_
- git push origin branchname _(push local commits/changes to the server)_
- git checkout master _(switch to master branch)_
- git pull _(pull down changes from master)_
- git checkout branchname, git merge master _(merge changes from master into your branch)_
- git status _(see information on the state of your local repository)_

If you have any questions about algorithms, architecture, implementations, or version control, feel free to email us at nwei@princeton.edu or kevinpg@princeton.edu!
