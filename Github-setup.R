

library(usethis)
use_git_config(user.name = "ondrejvladyka", user.email = "ondrejvladyka1@gmail.com")

usethis::use_git()

usethis::use_github()


gitcreds::gitcreds_set()

usethis::git_sitrep()


system('git config --global user.email "ondrejvladyka1@gmail.com"')
system('git config --global user.name "ondrejvladyka"')


system('git add Github-setup.R')

# 2. "Amend" the previous commit (this replaces the bad one with a clean version)
system('git commit --amend --no-edit')

# 3. Try to push again
system('git push')



system('git reset --soft origin/master')


system('git add .')
system('git commit -m "Cleaned scripts and added analysis"')
system('git push')



system('git add .')







system('git rm -r --cached .')









