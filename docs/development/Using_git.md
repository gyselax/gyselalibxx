# Using Git

Git is an extremely important tool for handling versioning and concurrent developments on the same code. There are many resources available online to help you get to grips with this tool. Here we summarise a few key concepts which will be helpful when working with git and Gyselalib++. It is designed to be a jumping-off point and a reference for reoccurring git questions, not a complete git reference or tutorial.

## Branches

Branches are an important tool for handling concurrent development. The `main` branch is the starting point for all developments. It contains all major developments which are considered sufficiently mature for other developers to also use them.

Whenever you wish to add something new to Gyselalib++ or fix an issue you should create a new branch from the `main` branch. The branch name should have the following format:

```none
<user_name>_<description_of_development>
```

this allows us to quickly and easily identify the branch owner and the contents of the branch.

It is important to keep your branch up to date with the main branch so that your changes can one day be merged into the main and used by other people. If this is not done on a regular basis then changes can accumulate making it exponentially harder to determine which changes are still relevant. As an **absolute minimum** you should update your branch **once a month**. You can do this using either a [merge](https://git-scm.com/docs/git-merge) command:

```sh
git fetch
git merge origin/main
```

or a [rebase](https://git-scm.com/docs/git-rebase) command:

```sh
git fetch
git rebase origin/main
```

After this command conflicts may appear between your branch and the main branch. For a `git merge` command all conflicts will appear at the same time. For a `git rebase` command conflicts may appear for each commit in your branch since the last rebase command (or its creation from the main branch). Conflicts should each be examined individually to ensure no new developments are lost. There are many resources (e.g. [bitbucket](https://www.atlassian.com/git/tutorials/using-branches/merge-conflicts)) online which describe how to fix these conflicts.

## Submodules

Gyselalib++ depends on other external libraries to compile. These libraries are included in the repository via [submodules](https://git-scm.com/docs/gitsubmodules). In general you will not need to touch the code inside the folders associated with these repositories. However their presence sometimes causes unexpected behaviour. The following is a small FAQ of common issues concerning submodules:

### Q: I cloned the repository but the submodules were not cloned

**A:** If you accidentally cloned the repository without the `--recurse-submodules` option you can collect the submodules by running the following command:

```sh
git submodule update --init
```

### Q: Git reports changes in the submodule but I didn't change this code

**A:** It is likely that you previously checked out a branch where the submodule was pinned to a different version. When you changed to your current branch the submodule was not updated, hence the reported changes. You can revert the changes to return to the version in your branch using:

```sh
git submodule update --init
```
