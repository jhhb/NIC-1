# NIC-1

#How to use git (at least for Mac)

#initializing
1. Create a new directory 'dir'
2. cd into dir
3. type "git init" - this initializes an empty github repository.
4. Keep this terminal window open

#cloning
1. Navigate to the github link for this portfolio on my page
2. Click the "Clone or download" green button and copy the link
3. Go to your terminal window and type "git clone copied_url" to clone the repository
4. You will now have a cloned repository here, meaning that you have the most recent version of the code at the time you clone.

#workflow
1. Before you start making changes, do 'git pull' in the directory to get the 
most recent code
2. Make your changes
3. Tell everyone you are pushing your changes and make sure this is ok with
everyone and that they dont have changes they need to push. If they do have changes they need to push, they should push them first.
4. Type 'git add --all'. This stages these changes for commit
5. Type 'git commit -m "descriptive message about changes goes here" '
6. Type 'git push origin master'
7. You repeat this workflow any time you make changes 

My instructions operate on the assumption that people won't be editing
The same part of the code. If people aren't editing the same code, it is easy to save your changes in a new text file window, update your local branch if someone else
needs to make changes first, and then put your changes into the updated code,
add, commit, and push

What we are doing, working on the master branch, is not great, and it
might end up being easier to use our own branches and merge.


