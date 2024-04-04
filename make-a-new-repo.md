## Gameplan for moving a GitHub repository subdirectory to its own new repository
### Relevant discussion with GPT4
<details>
<summary><i>Relevant discussion with GPT4</i></summary>
<br />

**Me**  
I want to move a GitHub repository subdirectory to its own new repository. Any advice regarding this?

**GPT4**  
Moving a subdirectory from one GitHub repository to its own new repository while preserving its commit history involves a few steps. Here's a step-by-step guide to help you through the process:

***Step 1: Clone the Original Repository***  
First, you'll need a local copy of the original repository. Replace `original-repo-url` with the URL of your GitHub repository.
```bash
git clone original-repo-url
cd name-of-original-repo
```

***Step 2: Filter the Subdirectory***  
Use the `git filter-branch` command or the newer `git-filter-repo` command (which you might need to install separately) to rewrite the repository's history, keeping only the commits that affect the specified subdirectory.

Using `git filter-branch` (pre-installed with Git but considered somewhat legacy):
```bash
git filter-branch --prune-empty --subdirectory-filter path/to/subdirectory HEAD
```

Using `git-filter-repo` (faster and more flexible, but requires separate installation):
First, install `git-filter-repo` if you haven't already (see the official documentation for installation instructions).

Then, run:
```bash
git-filter-repo --subdirectory-filter path/to/subdirectory
```

This command rewrites the repository's history to make the root of the subdirectory the new root of the repository.

***Step 3: Create a New Repository on GitHub***  
Go to GitHub and create a new repository for the subdirectory. Note the URL of the new repository.

***Step 4: Change the Remote Repository URL***  
Now, you'll want to push the filtered local repository to the new GitHub repository you just created. First, change the remote URL to the new repository's URL:
```bash
git remote set-url origin new-repo-url
```

Replace `new-repo-url` with the URL of the new repository you created on GitHub.

***Step 5: Push to the New Repository***  
Finally, push the local repository's content to the new GitHub repository:
```bash
git push -u origin master
```

If you're using a branch name other than `master` (such as `main` or something else), adjust the command accordingly.

***Additional Considerations***  
- *Large Repositories:* For very large repositories, the filter-branch command can be slow. In such cases, git filter-repo is recommended due to its better performance.
- *Dependencies and References:* If the subdirectory has dependencies on other parts of the original repository or if there are references you need to update (such as CI configurations), make sure to adjust those after moving the subdirectory.
- *Archiving:* If the original subdirectory will no longer be used in the original repository, consider removing it from there and committing that change to avoid confusion.

By following these steps, you can effectively move a subdirectory from one GitHub repository to its own new repository while retaining its commit history.
</details>
<br />
<br />

## On installing `git-filter-repo`
### Relevant conversation with GPT4
<details>
<summary><i>Relevant discussion with GPT4</i></summary>
<br />

**Me**  
To install `git-filter-repo`, I just need to put the script somewhere in my system (e.g., maybe my `${HOME}`) and then add it to my `${PATH}`. How should I go about this?

**GPT4**  
To install `git-filter-repo` and make it accessible from anywhere on your system, you indeed need to place the script in a directory that's included in your system's `${PATH}` variable. Here's a step-by-step guide on how to do this on a Unix-like system (Linux/macOS):

***Step 1: Download `git-filter-repo`***  
First, download `git-filter-repo` from its official repository. You can do this by visiting the `git-filter-repo` GitHub page and downloading the script manually or by using `wget` or `curl` in the terminal. For example:

```bash
curl -O https://raw.githubusercontent.com/newren/git-filter-repo/main/git-filter-repo
```

***Step 2: Move the Script to a Suitable Location***  
It's a good practice to keep custom scripts in a dedicated directory within your home directory, such as `${HOME}/bin`. If this directory doesn't exist, create it:
```bash
mkdir -p ${HOME}/bin
```

Then move `git-filter-repo` to this directory:
```bash
mv git-filter-repo ${HOME}/bin/
```

***Step 3: Make the Script Executable***  
Change the permissions of the script to make it executable:
```bash
chmod +x ${HOME}/bin/git-filter-repo
```

***Step 4: Add the Directory to Your `${PATH}`***  
To make sure your system can find `git-filter-repo` when you run it from any location, add `${HOME}/bin` to your `${PATH}`. Edit your shell's configuration file (`.bashrc` for Bash, `.zshrc` for Zsh, etc.):

For Bash, edit `.bashrc`:
```bash
nano ${HOME}/.bashrc
```

For Zsh, edit `.zshrc`:
```bash
nano ${HOME}/.zshrc
```

Add the following line to the end of the file:
```bash
export PATH="${HOME}/bin:${PATH}"
```

Save and close the file. Then, apply the changes:

```bash
source ${HOME}/.bashrc
```

Or, if you're using Zsh:
```bash
source ${HOME}/.zshrc
```

***Step 5: Verify the Installation***  
Verify that `git-filter-repo` is correctly installed and accessible from anywhere by running:

```bash
git-filter-repo --version
```

You should see the version number of `git-filter-repo` printed to the console.

This completes the installation. You now have `git-filter-repo` installed and added to your `${PATH}`, making it accessible from any location in your terminal.
</details>
<br />

## On using `usr/local/bin` as the location to store `git-filter-rep`
### Breakdown of what's in my `${PATH}`
Here's a breakdown of what's in my `${PATH}` as of 2024-0311:
<details>
<summary><i>Breakdown of what's in my `${PATH}`</i></summary>

```
/opt/homebrew/opt/grep/libexec/gnubin
/opt/homebrew/opt/gnu-which/libexec/gnubin
/opt/homebrew/opt/gnu-tar/libexec/gnubin
/opt/homebrew/opt/gnu-sed/libexec/gnubin
/opt/homebrew/opt/gnu-indent/libexec/gnubin
/opt/homebrew/opt/findutils/libexec/gnubin
/opt/homebrew/opt/ed/libexec/gnubin
/opt/homebrew/opt/libressl/bin
/opt/homebrew/opt/bison/bin
/opt/homebrew/opt/flex/bin
/opt/homebrew/opt/python/libexec/bin
/opt/homebrew/opt/unzip/bin
/opt/homebrew/opt/file-formula/bin
/opt/homebrew/opt/m4/bin
/opt/homebrew/opt/make/libexec/gnubin
/opt/homebrew/opt/coreutils/libexec/gnubin
/opt/homebrew/bin
/Users/kalavattam/miniconda3/bin
/Users/kalavattam/miniconda3/condabin
/opt/homebrew/bin
/opt/homebrew/sbin
/usr/local/bin
/System/Cryptexes/App/usr/bin
/usr/bin
/bin
/usr/sbin
/sbin
/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/local/bin
/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/bin
/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/appleinternal/bin
/opt/X11/bin
/opt/homebrew/opt/fzf/bin
/Users/kalavattam/GTFtools_0.8.0
/Users/kalavattam/bbmap
```
</details>
<br />

### Relevant discussion with GPT4
<details>
<summary><i>Relevant discussion with GPT4</i></summary>
<br />

**Me**  
Here are the `bin` options I have in my `$PATH`:  
`/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin`

Which one do you think I should work with for `git-filter-repo`?

**GPT4**  
Given the options in your `$PATH`, `/usr/local/bin` is generally the best choice for installing user-level software and scripts like `git-filter-repo`. This directory is intended for software and scripts that are not managed by the system's package manager and are intended to be accessible to all users on the system. Installing custom or third-party software in `/usr/local/bin` helps keep them separate from the system's default binaries located in `/usr/bin` and `/bin`, and it's a common practice on Unix-like systems.

Here's why `/usr/local/bin` is a suitable choice:
- *User-installed Software:* It's a conventional location for software and scripts installed by the user, rather than software provided by the system's package manager.
- *System-wide Accessibility:* Placing the script in /usr/local/bin makes it accessible to all users on the system (assuming they all have `/usr/local/bin` in their `$PATH`), which is useful for tools like `git-filter-repo` that you might want available regardless of the active user.
- *Avoids Overwriting System Binaries:* By using `/usr/local/bin`, you reduce the risk of overwriting system binaries or creating conflicts with system-managed files in `/usr/bin` and `/bin`.

***Steps to install `git-filter-repo` in `/usr/local/bin`***  
*Download `git-filter-repo` (if you haven't already):*
```bash
curl -O https://raw.githubusercontent.com/newren/git-filter-repo/main/git-filter-repo
```
*Move and rename (optional) the script to `/usr/local/bin`:*
```bash
sudo mv git-filter-repo /usr/local/bin/git-filter-repo
```
*Make it executable:*
```bash
sudo chmod +x /usr/local/bin/git-filter-repo
```
*Verify installation:*
```bash
git-filter-repo --version
```

By following these steps, you'll install `git-filter-repo` in a standard, system-wide accessible location that aligns with Unix conventions for user-installed software.
</details>
<br />
<br />

## What is the `new-repo-url` for a repository I have created?
### Relevant discussion with GPT4
<details>
<summary><i>Relevant discussion with GPT4</i></summary>
<br />

**Me**  
How can I figure out what is the `new-repo-url` for a repository I have created?

E.g., with respect to this command:
```bash
git remote set-url origin new-repo-url
```

**GPT4**  
To figure out the `new-repo-url` for a repository you've created, you'll need to obtain the URL provided by GitHub (or another Git hosting service) for that repository. Here's how you can do it on GitHub, which is a common scenario:

***Finding the Repository URL on GitHub***  
1. *Navigate to Your GitHub Repository:* Open your web browser and go to your GitHub account. Then, navigate to the repository you've created.
2. *Locate the "Code" Button:* On the repository's main page, you'll see a green button labeled "Code" near the top right, just above the file list.
3. *Get the Repository URL:* Click the "Code" button. A dropdown menu will appear, offering you options to clone the repository with HTTPS, SSH, or GitHub CLI. Here are your options:
- *HTTPS:* If you're using HTTPS, simply copy the URL shown (it will look something like https://github.com/username/repo-name.git).
- *SSH:* If you're using SSH, switch to the SSH tab and copy the SSH URL (it will look like git@github.com:username/repo-name.git).
- *GitHub CLI:* For using the GitHub CLI, it shows a command rather than a URL, so you'd typically use HTTPS or SSH URLs for setting remote URLs.

***Using the Repository URL in Git***  
Once you've copied the repository URL, use it in your command to set the new remote URL:
```bash
git remote set-url origin <new-repo-url>
```

Replace <new-repo-url> with the URL you just copied. For example:
```bash
git remote set-url origin https://github.com/username/repo-name.git
```

or with SSH:
```bash
git remote set-url origin git@github.com:username/repo-name.git
```

***Confirming the Change***  
After changing the remote URL, you can confirm that the operation was successful by listing the remote URLs associated with your repository:
```bash
git remote -v
```

This command will show you the URLs for fetching and pushing for each remote associated with your local repository. After executing the `set-url` command, you should see the remote URL for origin has been updated to the new URL you provided.

This process applies to GitHub, but the steps are similar for other Git hosting services like GitLab or Bitbucket, with differences mainly in how they structure their repository URLs and where the clone URL can be found on their respective web interfaces.
</details>
<br />
<br />

## On where to run `git remote set-url origin new-repo-url`
### Relevant discussion with GPT4
<details>
<summary><i>Relevant discussion with GPT4</i></summary>
<br />

**Me**  
Currently, I am in the following directory,
```txt
/Users/kalavattam/repos/bak/2023_rDNA
```

and I have just run the following command,
```bash
git-filter-repo --subdirectory-filter "results/2023-0406_tutorial_ChIP-seq_analyses"
```

which printed the following to my console,
```txt
Parsed 229 commits
New history written in 0.48 seconds; now repacking/cleaning...
Repacking your repo and cleaning out old unneeded objects
HEAD is now at ecf547c Finally pushing the last couple weeks of work (but not the scratch work) to the tutorial repo; this work concerned various set operations with ChIP-seq MACS3 peak sets (Rachel's data: Esa1, Gcn5, and Rpd3)
Enumerating objects: 170, done.
Counting objects: 100% (170/170), done.
Delta compression using up to 8 threads
Compressing objects: 100% (98/98), done.
Writing objects: 100% (170/170), done.
Total 170 (delta 76), reused 130 (delta 72), pack-reused 0
Completely finished after 0.70 seconds.
```

According to your instructions, the next step is as follows:
> Step 4: Change the Remote Repository URL  
> Now, you'll want to push the filtered local repository to the new GitHub repository you just created. First, change the remote URL to the new repository's URL:
> ```bash
> git remote set-url origin new-repo-url
> ```
> 
> Replace `new-repo-url` with the URL of the new repository you created on GitHub.

Am I to run `git remote set-url origin new-repo-url`, where `new-repo-url` is `https://github.com/kalavattam/2023_teaching.git` from my current directory, `/Users/kalavattam/repos/bak/2023_rDNA`?
</details>
<br />
<br />

## Install `git-filter-rep` to `/usr/local/bin`
### Relevant discussion with GPT4
<details>
<summary><i>Relevant discussion with GPT4</i></summary>
<br />

**Me**  
```bash
cd /usr/local/bin

sudo curl -O https://raw.githubusercontent.com/newren/git-filter-repo/main/git-filter-repo

sudo chmod +x git-filter-repo

git-filter-repo --version
# ae71fad1d03f
```

**GPT4**  
Yes, based on your current directory (`/Users/kalavattam/repos/bak/2023_rDNA`) and the steps you've followed so far with `git-filter-repo`, you should run the `git remote set-url origin new-repo-url` command with `new-repo-url` being https://github.com/kalavattam/2023_teaching.git from your current directory. This command will change the remote URL of the repository located in your current directory to point to the new repository on GitHub, allowing you to push the filtered history to the new repository.

Here's how to do it:

1. Ensure you are in the correct directory (`/Users/kalavattam/repos/bak/2023_rDNA`). If you're not sure, you can verify your current directory by running:
```bash
pwd
```

Change the remote URL to the new repository's URL:
```bash
git remote set-url origin https://github.com/kalavattam/2023_teaching.git
``` 

You can confirm that the remote URL was successfully updated by listing the current remotes with their URLs:
```bash
git remote -v
```
You should see the URL for the origin remote has been updated to https://github.com/kalavattam/2023_teaching.git.

After updating the remote URL, you're ready to push the changes to the new repository. Since `git-filter-repo` has rewritten the history, it is a good practice to use `git push --force` to overwrite the history in the remote repository with your newly filtered history:
```bash
git push --force origin master
```

If you're using a different default branch (such as `main`), replace `master` with the name of your default branch in the command above.
</details>
<br />
<br />

## What if there is no such remote as `origin`?
### Relevant discussion with GPT4
<details>
<summary><i>Relevant discussion with GPT4</i></summary>
<br />

**Me**  
How can I troubleshoot this? How can I find the name I should use instead of `origin`?
```txt
❯ git remote set-url origin https://github.com/kalavattam/2023_teaching.git
error: No such remote 'origin'
```

**GPT4**  
The error message `No such remote 'origin'` indicates that there is no remote repository with the name `origin` currently associated with your local Git repository. You can troubleshoot and resolve this issue by following these steps:

***List Existing Remotes***  
First, check the current remotes associated with your repository to see what they're named:
```bash
git remote -v
```

This command will list all the remote repositories your local repository knows about, along with the URLs associated with the fetch and push actions for each remote.

***Add a New Remote***  
If the output of `git remote -v` does not list an `origin` (or if there are no remotes listed at all), you need to add a new remote with the name `origin` pointing to your GitHub repository URL:
```bash
git remote add origin https://github.com/kalavattam/2023_teaching.git
```

***Rename an Existing Remote***  
If there is a remote, but it's not named `origin` (for example, if it's named upstream or some other name), and you want to rename it to `origin`, you can use the `git remote rename` command:
```bash
git remote rename existing_remote_name origin
```

Replace `existing_remote_name` with the actual name of the remote you found using `git remote -v`.

***Set the URL for an Existing Remote***  
If there's already a remote with a different name, and you simply want to change its URL without renaming it to `origin`, you can replace `origin` in the original command with the actual name of your remote:
```bash
git remote set-url existing_remote_name https://github.com/kalavattam/2023_teaching.git
```

Again, replace `existing_remote_name` with the name of the remote you're intending to update.

***Confirm the Changes***  
After adding, renaming, or updating the remote, run `git remote -v` again to confirm that the changes were made successfully. You should now see the origin remote (or whichever name you've chosen to use) listed with the correct URL to your GitHub repository.

With the remote correctly set up, you can proceed with pushing your changes to the repository using `git push` (you might need to use `--force` or `--set-upstream` options depending on your situation).
</details>
<br />
<br />

## Is it enough to run `git remote add origin` prior to running `git push -u origin master`?
### Relevant discussion with GPT4
<details>
<summary><i>Relevant discussion with GPT4</i></summary>
<br />

**Me**  
I have run the following command:
```bash
git remote add origin https://github.com/kalavattam/2023_teaching.git
```

Having done so, can I proceed with Step #5?
> Step 5: Push to the New Repository
> Finally, push the local repository's content to the new GitHub repository:
> 
> ```bash
> git push -u origin master
> ```
> If you're using a branch name other than master (such as main or something else), adjust the command accordingly.

**GPT4**  
Yes, after successfully adding a remote named `origin` with the command you've run:
```bash
git remote add origin https://github.com/kalavattam/2023_teaching.git
```

You can proceed with Step #5 and push the local repository's content to the new GitHub repository. However, you should check which branch you are currently on, especially since many repositories now use `main` as the default branch name instead of `master`.

To check your current branch, use:
```bash
git branch
```

The branch you're currently on will be marked with an asterisk \(\*\).

If your current branch is `master`, you can follow the step as written:
```bash
git push -u origin master
```

If your branch is named differently, such as main or any other name, you should replace master with the name of your current branch:
```bash
git push -u origin main
```

The `-u` flag sets the upstream for your branch, linking your local branch with the remote one, which will make future pushes (simply `git push`) and pulls automatically refer to this upstream setting without needing to specify the remote and branch names.

If you have just filtered your repository and it's the first time you are pushing to this new repository, you may need to use the `--force` flag because `git-filter-repo` rewrites the history, which can cause discrepancies with the remote history:
```bash
git push -u origin main --force
```

Replace `main` with your actual branch name if it's different.
</details>
<br />
<br />

## On encountering an error with `git push -u origin main`
### Relevant discussion with GPT4
<details>
<summary><i>Relevant discussion with GPT4</i></summary>
<br />

**Me**  
After running this command,
```bash
git push -u origin main
```

I got the following error message:
```txt
❯ git push -u origin main
To https://github.com/kalavattam/2023_teaching.git
 ! [rejected]        main -> main (fetch first)
error: failed to push some refs to 'https://github.com/kalavattam/2023_teaching.git'
hint: Updates were rejected because the remote contains work that you do
hint: not have locally. This is usually caused by another repository pushing
hint: to the same ref. You may want to first integrate the remote changes
hint: (e.g., 'git pull ...') before pushing again.
hint: See the 'Note about fast-forwards' in 'git push --help' for details.
```

Regarding `2023_teaching`, I don't care about any "work that \[I\] do not have locally". How can I move forward? Is it possible for me to force the `push`?

**GPT4**  
Yes, you can force the push if you're sure that you don't need any commits that are currently on the remote repository and that overwriting the history on the remote is safe in your context. This situation often arises after using tools like `git-filter-repo`, where you've intentionally rewritten the history of your repository and need the remote repository to reflect these changes.

To force the push to the remote repository and overwrite its history with your local history, you can use the `--force` flag with `git push`:
```bash
git push -u origin main --force
```

Or, as a safer alternative, you can use `--force-with-lease`, which ensures you do not overwrite the work of others on the remote repository in case there have been updates to it since you last fetched:
```bash
git push -u origin main --force-with-lease
```

The `--force-with-lease` option is generally recommended over `--force` because it provides an extra level of safety, ensuring that you do not accidentally overwrite remote changes you're unaware of. However, if you're absolutely certain that it's safe to overwrite the remote repository (for instance, if it's a new repository or you're the only one working on it), using `--force` is acceptable.

Remember, using `--force` or `--force-with-lease` will overwrite the history on the remote. This action cannot be undone and may result in loss of work if not used carefully. Always make sure that this is what you really intend to do before proceeding.
</details>
<br />
<br />

## On actually filtering the repository subdirectory&mdash;i.e., moving a GitHub repository subdirectory to its own new repository
### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

cd /usr/local/bin

sudo curl -O https://raw.githubusercontent.com/newren/git-filter-repo/main/git-filter-repo

sudo chmod +x git-filter-repo

git-filter-repo --version

cd ~/repos

mkdir bak

cd bak

gh repo clone kalavattam/2023_rDNA

cd 2023_rDNA

git-filter-repo --subdirectory-filter "results/2023-0406_tutorial_ChIP-seq_analyses"

git remote set-url origin https://github.com/kalavattam/2023_teaching.git
#  Error

git remote -v
#  Returns nothing

git remote add origin https://github.com/kalavattam/2023_teaching.git
#  Seems to have worked 

git branch
# * main

# git push -u origin main  # See question below
git push -u origin main --force
# Enumerating objects: 170, done.
# Counting objects: 100% (170/170), done.
# Delta compression using up to 8 threads
# Compressing objects: 100% (94/94), done.
# Writing objects: 100% (170/170), 2.97 MiB | 3.40 MiB/s, done.
# Total 170 (delta 76), reused 170 (delta 76), pack-reused 0
# remote: Resolving deltas: 100% (76/76), done.
# To https://github.com/kalavattam/2023_teaching.git
#  + 6236e05...ecf547c main -> main (forced update)
# branch 'main' set up to track 'origin/main'.
```
</details>
<br />
