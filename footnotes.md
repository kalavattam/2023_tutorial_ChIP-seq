
`#footnotes.md`
<br />

<details>
<summary><b><font size="+2"><i>Table of contents</i></font></b></summary>
<br />
<!-- MarkdownTOC -->

1. [1. Prepare FASTQ files of interest for processing and analyses](#1-prepare-fastq-files-of-interest-for-processing-and-analyses)
	1. [Footnotes](#footnotes)
		1. [Comments \(`#`\)](#comments-)
		1. [Defining a function such as `error_and_return`](#defining-a-function-such-as-error_and_return)
		1. [The body of the function `error_and_return`](#the-body-of-the-function-error_and_return)
		1. [Variables](#variables)
		1. [The `${HOME}` variable](#the-%24home-variable)
		1. [Wrapping variables in curly braces `{}` and double quotes `""`](#wrapping-variables-in-curly-braces--and-double-quotes-)
		1. [Associative arrays \(hash maps\)](#associative-arrays-hash-maps)
		1. [Logical commands \(`true`, `false`\) assigned to "flag variables" \(flags\)](#logical-commands-true-false-assigned-to-flag-variables-flags)
		1. [`if` statements](#if-statements)
		1. [Conditional checks for directories \(`-d`\), files \(`-f`\), and logical negation \(`!`\)](#conditional-checks-for-directories--d-files--f-and-logical-negation-)
		1. [Calls to `ln`](#calls-to-ln)
		1. [Calls to `ls`](#calls-to-ls)
1. [2. Adapter- and quality-trim the FASTQ files](#2-adapter--and-quality-trim-the-fastq-files)
	1. [a. Install Atria and dependencies](#a-install-atria-and-dependencies)
		1. [Footnotes](#footnotes-1)
			1. [Breaking down the function `update_shell_config`](#breaking-down-the-function-update_shell_config)
				1. [Positional arguments](#positional-arguments)
				1. [Local variable scoping](#local-variable-scoping)
				1. [Redirection: `>>` versus `>`](#redirection--versus-)
				1. [Calls to `grep`](#calls-to-grep)
				1. [Return values](#return-values)
			1. [Breaking down the function `check_mamba_installed`](#breaking-down-the-function-check_mamba_installed)
				1. [`if` statement with negation](#if-statement-with-negation)
				1. [Redirection: `&> /dev/null`](#redirection--devnull)
				1. [Redirection: `&>`, `1>`, `2>`, and more](#redirection--1-2-and-more)
				1. [On data streams: stdin, stdout, stderr](#on-data-streams-stdin-stdout-stderr)
			1. [Breaking down the function `check_env_installed`](#breaking-down-the-function-check_env_installed)
				1. [Declaration of the function](#declaration-of-the-function)
				1. [Local variable declaration, e.g., `local env_name="${1}"`](#local-variable-declaration-eg-local-env_name%241)
				1. [The Conda/Mamba command `conda env list`](#the-condamamba-command-conda-env-list)
				1. [The pipe \(`|`\) and `grep` commands: `... | grep -q "^${env_name} "`](#the-pipe-%7C-and-grep-commands--%7C-grep--q-%5E%24env_name-)
				1. [Regular expressions and the caret \(`^`\) symbol:](#regular-expressions-and-the-caret-%5E-symbol)
				1. [Escape characters and `\"`:](#escape-characters-and-)
				1. [On the function's "control flow"](#on-the-functions-control-flow)
	1. [b. Adapter- and quality-trim the FASTQ files using Atria](#b-adapter--and-quality-trim-the-fastq-files-using-atria)
		1. [Footnotes](#footnotes-2)
			1. [Calls to `unset`](#calls-to-unset)
			1. [Calls to `typeset`/`declare`](#calls-to-typesetdeclare)
			1. [Common options for `typeset`/`declare`](#common-options-for-typesetdeclare)
			1. [More on operators, particularly the logical operators `&&` and `||`](#more-on-operators-particularly-the-logical-operators--and-%7C%7C)
				1. [`&&` Operator \(Logical AND\)](#-operator-logical-and)
				1. [`||` Operator \(Logical OR\)](#%7C%7C-operator-logical-or)
				1. [How `&&` and `||` work together](#how--and-%7C%7C-work-together)
				1. [Other relevant operators: `;` and `!` \(logical NOT\)](#other-relevant-operators--and--logical-not)
				1. [Best Practices](#best-practices)
				1. [Summary](#summary)
1. [3. Align trimmed FASTQ files](#3-align-trimmed-fastq-files)
	1. [a. Generate a concatenated annotated assembly of the *S. cerevisiae* and *S. pombe* genomes](#a-generate-a-concatenated-annotated-assembly-of-the-s-cerevisiae-and-s-pombe-genomes)
		1. [Footnotes](#footnotes-3)
			1. [Breaking down the call to `mkdir -p`, which makes use of brace expansion](#breaking-down-the-call-to-mkdir--p-which-makes-use-of-brace-expansion)

<!-- /MarkdownTOC -->
</details>
<br />

<a id="1-prepare-fastq-files-of-interest-for-processing-and-analyses"></a>
# 1. Prepare FASTQ files of interest for processing and analyses
<a id="footnotes"></a>
## Footnotes
<details>
<summary><i>Notes: 1. Prepare FASTQ files of interest for processing and analyses</i></summary>
<br />

<a id="comments-"></a>
### Comments (`#`)
The lines starting with `#` are comments. Comments are used to explain what the code does, but they are not executed. They're there to help anyone reading the code understand it better.

For example, `# Define functions ===========================================================` is a header comment, letting readers know that functions are defined below it. The next comment explains what the function `error_and_return` does.

<a id="defining-a-function-such-as-error_and_return"></a>
### Defining a function such as `error_and_return`
`function error_and_return()` starts the definition of a function named `error_and_return`. A function is a block of code that performs a specific task. Functions facilitate the reuse of the same code multiple times without having to write it out each time.

<a id="the-body-of-the-function-error_and_return"></a>
### The body of the function `error_and_return`
Inside the function, we have the following code:
- `echo "Error: ${1}" >&2`
    + `echo` is a command used to display text.
    + `"Error: ${1}"` is the text to be displayed. Here, `${1}` is a placeholder for the first argument passed to the function. This means that whatever text is provided when calling `error_and_return` will be displayed after `"Error: "`.
    + `>&2` means that the error message is redirected to the standard error stream (stderr). This is typically used to output error messages.
- `return 1`
    + `return 1` exits the function and returns a value of 1. In Unix, Linux, and MacOS, returning a non-zero value generally indicates an error or abnormal condition. So, when `error_and_return` is called, it indicates that something went wrong.

<a id="variables"></a>
### Variables
Variables (e.g., `dir_base`) are used to store and retrieve data (e.g., `dir_base="${HOME}/tsukiyamalab"`). They are like placeholders for values that can change over time. Variables make scripts flexible and reusable. For example, changing `dir_base` (to, say, `dir_base="${HOME}/Desktop"`) will update the base directory used throughout the script.

<a id="the-%24home-variable"></a>
### The `${HOME}` variable
The `${HOME}` variable is an environmental variable that refers to the home directory of the current user. In Unix and Unix-like operating systems (e.g., MacOS), every user is assigned a unique directory where they can store personal files, configurations, and scripts. This directory is commonly referred to as the "home directory."
- In scripts and command lines, `${HOME}` is often used as a shorthand to access the user's home directory. For example, `cd ${HOME}` would change the current directory to the user's home directory.
- If the username is `john`, the home directory might be `/home/john` on Linux or `/Users/john` on macOS. In this case, `${HOME}` would be equivalent to `/home/john` or `/Users/john`, respectively.
- Using `${HOME}` makes scripts more portable and user-independent, as it automatically resolves to the home directory of the user running the script, without needing to hard-code the full path.

<a id="wrapping-variables-in-curly-braces--and-double-quotes-"></a>
### Wrapping variables in curly braces `{}` and double quotes `""`
In shell scripting, it's considered good practice to wrap variables in curly braces (`{}`) and double quotes (`""`) for clarity and to prevent potential errors. This practice has several benefits:
- Curly braces help in clearly defining the boundary of a variable name. This is particularly useful when a variable is followed by text that could be misinterpreted as part of the variable name. For example, `${var}_suffix` clearly separates `var` from `_suffix`.
- Curly braces facilitate the concatenation of variables with strings. For instance, `${var}value` appends `value` to the contents of `var`.
- Enclosing variables in double quotes prevents word splitting and globbing. Word splitting can lead to unexpected behavior when variables contain spaces or newlines. For example, in `echo ${var}`, if `var` contains `file1 file2`, it will be split into two arguments. `echo "${var}"` would treat it as a single argument, preserving the intended behavior.
- When a variable is empty or undefined, using double quotes ensures that the script doesn't break or behave unpredictably. For example, `echo "${nonexistent}"` will safely print nothing, whereas `echo ${nonexistent}` might lead to unintended script behavior.
- Example:
    + With braces: `echo "${user}file"` makes it clear that `user` is the variable being referenced, and the string `"file"` is being appended to its contents.
    + Without braces: `echo "$userfile"` is misread as if there were a variable named `userfile`.

The use of curly braces and double quotes enhances the readability and reliability of shell scripts, making them less prone to errors.


<a id="associative-arrays-hash-maps"></a>
### Associative arrays (hash maps)
Associative arrays (or hash maps) are collections of key-value pairs where each key is unique. They allow for more complex data structures, enabling you to map a unique key to a specific value. In the above chunk, keys are original file name stems and values are the new name stems for the symbolic links. For example, `file_fastqs["6336_G1_in_S15"]="in_G1_Hho1_6336"` maps an original file name stem, `"6336_G1_in_S15"`, to a new one, `"in_G1_Hho1_6336"`.

<a id="logical-commands-true-false-assigned-to-flag-variables-flags"></a>
### Logical commands (`true`, `false`) assigned to "flag variables" (flags)
Flags are variables used to control the flow of the script. `check_variables`, `check_array`, `check_operations` and `run_operations` are flags. When set to `true`, they trigger specific operations like echoing commands or creating symbolic links. Conversely, setting them to `false` skips these operations.

<a id="if-statements"></a>
### `if` statements
`if` statements are used for the conditional execution of code. They allow the script to make decisions and execute different blocks of code based on whether a condition is true or false. They follow the following logic: "`if` a given condition is `true`, execute a specific code block; `else` (optionally) execute a different code block or do nothing." For example, `if [[ -d "${directory}" ]]; then echo "Directory exists"; fi`.

<a id="conditional-checks-for-directories--d-files--f-and-logical-negation-"></a>
### Conditional checks for directories (`-d`), files (`-f`), and logical negation (`!`)
- `[[ -d ... ]]`: Checks if a directory exists.
- `[[ -f ... ]]`: Checks if a file exists.
- `!`: Negates a condition, e.g., `[[ ! -d ... ]]` checks if a directory *does not* exist.

<a id="calls-to-ln"></a>
### Calls to `ln`
The `ln -s` command creates symbolic links (symlinks), which are pointers to original files. The command format is `ln -s original_file symlink_file`. Using `ln` like this maintains the integrity of raw data while simplifying access under new names.

<a id="calls-to-ls"></a>
### Calls to `ls`
The `ls` command lists directory contents. Flags used include:
- `-l`: Long format with detailed information.
- `-h`: Human-readable file sizes.
- `-a`: Includes hidden files.
- `-F`: Appends a character indicating file type.
- `-G`: Colorizes output.

For more details on `ls` flags and shell commands in general, visit [ShellCheck](https://www.shellcheck.net/).
</details>
<br />
<br />

<a id="2-adapter--and-quality-trim-the-fastq-files"></a>
# 2. Adapter- and quality-trim the FASTQ files
<a id="a-install-atria-and-dependencies"></a>
## a. Install Atria and dependencies
<a id="footnotes-1"></a>
### Footnotes
<details>
<summary><i>Notes: 2.a. Install Atria and dependencies</i></summary>
<br />

`#TODO` Carefully explain all concepts in the above chunk.  
`#TODO` Carefully test the code in the above chunk.

<a id="breaking-down-the-function-update_shell_config"></a>
#### Breaking down the function `update_shell_config`
<a id="positional-arguments"></a>
##### Positional arguments
When a function is called, you can pass data to it via "arguments". In this function, `${1}` is a placeholder for the first argument, i.e., the argument in "the first position"; `${2}` is a placeholder for the second argument passed, i.e., the argument in "the second position". The values passed to positional arguments `${1}` and `${2}` are assigned to the "local" variables `config_file` and `stem`, respectively.

<a id="local-variable-scoping"></a>
##### Local variable scoping
The `local` command is used to declare variables that are "local" to the function. This means these variables (`config_file`, `stem`, `line_to_add`) only exist within `update_shell_config` and can't be accessed outside of it. This means the "scope" of these variables are "local" to the function. Local variable scoping helps prevent conflicts with variables of the same name elsewhere in the script.

<a id="redirection--versus-"></a>
##### Redirection: `>>` versus `>`
- `>>` appends the output (stdout) of a command or operation to a file. If the file doesn't exist, then it's created. If it does exist, then the new output is added at the end of the file.
- Contrast that with `>`, which redirects output to a file, overwriting its current contents. Again, if the file doesn't exist, then it's created.

<a id="calls-to-grep"></a>
##### Calls to `grep`
`grep` is a command-line utility for searching plain-text data for lines that match a regular expression. Here, `grep -q "${line_to_add}" "${config_file}"` checks if the `line_to_add` is already in `config_file`. The `-q` flag makes `grep` operate in quiet mode, so it doesn't output anything and instead just returns a success (`0`) or failure (`1`) status.

<a id="return-values"></a>
##### Return values
In shell scripting, the command `return` exits the function. `return 0` typically signifies success, and `return 1` (or any non-zero value) signifies failure or an error. These values can be used by other parts of the script to determine if the function succeeded or failed.

<a id="breaking-down-the-function-check_mamba_installed"></a>
#### Breaking down the function `check_mamba_installed`
This function does the following: `#TODO`.

<a id="if-statement-with-negation"></a>
##### `if` statement with negation
- `if` statements (e.g., see the code starting with `if ! type mamba &> /dev/null; then`) check a condition and execute code based on whether the condition is true or false.
- The `!` before a command negates its success status.
- So, if `type mamba &> /dev/null` is successful and returns an exit code of 0, then the `!` in `! type mamba &> /dev/null` would convert that exit code of `0` to an exit code of `1`.
- With that in mind, checking for the mamba command works in the following way:
    + `! type mamba &> /dev/null` returns `0` if mamba is not actually in the PATH (as `!` converts `1` to `0`).
    + Otherwise, `! type mamba &> /dev/null` returns `1` (as `!` converts `0` to `1`), thereby skipping the block of code within the `if` statement and leading to the `return 0` command.

<a id="redirection--devnull"></a>
##### Redirection: `&> /dev/null`
`&> /dev/null` redirects both the standard output (stdout) and standard error (stderr) to `/dev/null`, effectively silencing all the command's output. `/dev/null` is a special file unique to Unix and Unix-like operating systems that discards all data written to it. This makes `/dev/null` useful for suppressing unwanted output from commands or scripts.

<a id="redirection--1-2-and-more"></a>
##### Redirection: `&>`, `1>`, `2>`, and more
In shell scripting, redirection operators are used to control where the output of commands goes. 
- Standard output redirection: `>` or `1>`
    + This operator redirects the standard output (stdout) of a command to a file or another command.
    + In shell scripting, `1` represents the file descriptor for standard output. Since it's the default, `>` is equivalent to `1>`.
    + For example, `echo "Hello, World" > file.txt` writes "Hello, World" to `file.txt`.
- Standard error redirection: `2>`
    + This operator redirects the standard error (stderr) of a command to a file or another command.
    + In shell scripting, `2` represents the file descriptor for standard error.
    + For example, `ls non_existent_file 2> error.txt` redirects any error messages from the `ls` command to `error.txt`.
- Combined stdout and stderr redirection: `&>`
    + This operator redirects both standard output and standard error to the same place.
    + It is useful to capture all output from a command, regardless of whether it is normal output or error messages.
    + For example, `command &> output.txt` will redirect both the output and any error messages of `command` to `output.txt`.
- Appending stdout redirection: `>>`
    + This operator appends the standard output of a command to the end of an existing file, rather than overwriting it like `>` does.
    + For example, `echo "World" >> file.txt` will add "World" to the end of `file.txt` without removing any existing content.
- Redirecting stderr to stdout: `2>&1`
    + This operator redirects the standard error to the same destination as the standard output.
    + `2>&1` is often used in combination with other redirections; for example, `command > output.txt 2>&1` will redirect both stdout and stderr to `output.txt`.

<a id="on-data-streams-stdin-stdout-stderr"></a>
##### On data streams: stdin, stdout, stderr
In computing, particularly in the context of Unix and Unix-like operating systems, data streams are channels through which data flows. The three standard streams are:
1. Standard input (stdin): This is the data stream used for input. Typically, stdin is what you type into the terminal. By default, stdin is "attached to" or associated with the keyboard. To visualize this, imagine a natural stream (body of water) whose source is underground water; similarly, the keyboard acts as the source of "signals" or "data" (like the water), initiating a "data stream" (stdin). This stream "flows" into programs or commands, carrying the input (data) they require.
2. Standard output (stdout): This stream is used to output the data produced by a program. For example, when you run a command in the terminal that prints something, that output is being sent to stdout. By default, stdout is displayed on the screen. In the context of a terminal, the river's banks are the screen where output of commands are viewed. In engineering, people might direct a river through channels or pipes to specific locations. Similarly, stdout can be redirected to files, other programs, or even other devices. This redirection is akin to building a canal or pipeline to guide the river's flow to a desired destination.
3. Standard error (stderr): This is a separate stream used specifically for outputting error messages or diagnostics from a program. It is distinct from stdout, which allows error messages to be handled or redirected separately from standard output. By default, stderr is also displayed on the screen. Engineers often design separate drainage systems to handle waste or overflow. Similarly, stderr can be redirected independently of stdout. This is like having a separate set of pipes or channels (e.g., a sewage system) to manage waste water, ensuring it doesn’t pollute the main water supply.
4. Taking it all in:
    + For example, in nature, water can be directed using various natural formations like canyons and gorges, or man-made structures like dams, sluice gates, and reservoirs. In "\*nix" operating systems (e.g., Unix, Linux, MacOS), operators like `>`, `>>`, `2>`, `|`, and others are used as tools to redirect and control these data streams. They act like, for example, sluice gates or pipes (`|`) (to redirect flow), dams (to stop and store data), or even filters (to process and change data).
    + So, the stdin, stdout, and stderr streams are fundamental in \*nix environments for data input and output. Just as water can be filtered, stored, or channeled into different paths (like irrigation systems, hydroelectric plants, or through filtration systems), data streams can be manipulated and redirected in numerous ways. The use of commands and scripts to manipulate these streams is akin to using sophisticated control systems in engineering to manage water flow, ensuring each drop goes exactly where it's needed, when it's needed.

<a id="breaking-down-the-function-check_env_installed"></a>
#### Breaking down the function `check_env_installed`
This function checks for the existence of a specified Conda environment. To achieve its goal, it uses "local variable scoping" (`local env_name="${1}"`) with a single positional argument (`${1}`), a conditional statement (`if`, `then`, `else`), pattern matching (`grep -q`), a regular expression (the `^` in `"^${env_name} "`), and command piping (`|`).

<a id="declaration-of-the-function"></a>
##### Declaration of the function
`function check_env_installed() { ... }` defines a new function named `check_env_installed`
<a id="local-variable-declaration-eg-local-env_name%241"></a>
##### Local variable declaration, e.g., `local env_name="${1}"`
- The `local` command restricts the variable's scope to the function.
- `env_name` is a variable that holds the name of the environment to check.
- `"${1}"` is the first positional argument passed to the function when it's called.

<a id="the-condamamba-command-conda-env-list"></a>
##### The Conda/Mamba command `conda env list`
`conda env list` command lists all Conda environments installed on the system.

<a id="the-pipe-%7C-and-grep-commands--%7C-grep--q-%5E%24env_name-"></a>
##### The pipe (`|`) and `grep` commands: `... | grep -q "^${env_name} "`
- The pipe (`|`) takes the output of `conda env list` and passes it to the `grep` command, which searches for a specific pattern in the input it receives.
- `-q` is an option for `grep` that makes it silent; with `-q` specified, `grep` doesn't output the matching lines; it returns an exit status (`0` or a non-zero value).
- `"^${env_name} "` is the pattern `grep` searches for.

<a id="regular-expressions-and-the-caret-%5E-symbol"></a>
##### Regular expressions and the caret (`^`) symbol:
- The pattern `"^${env_name} "` includes a regular expression. Regular expressions are a way to match patterns in text.
- `^` is a regular expression "anchor" that matches the start of a line.
- `^${env_name}` means `grep` looks for lines that start with the name of the environment&mdash;the value assigned to variable `env_name`.

<a id="escape-characters-and-"></a>
##### Escape characters and `\"`:
- The backslash (`\`) is used as an "escape character". It changes the meaning of the character following it.
- Here, it is used to include the double quotes literally in the output string, i.e., `echo "Environment \"${env_name}\" is not installed."` prints the environment name within quotes.

<a id="on-the-functions-control-flow"></a>
##### On the function's "control flow"
- If the specified environment is found (that is, if the pattern is matched), `grep -q` returns a zero exit status (`0`), which indicates success. Consequently, the function also returns `0`.
- If the environment is not found (that is, if the pattern is not matched), `grep -q` returns a non-zero status. This means that the `else` block (i.e., the block of code following (below) the `else` operator) executes, printing a message and returning `1`, which indicates failure.
</details>
<br />

<a id="b-adapter--and-quality-trim-the-fastq-files-using-atria"></a>
## b. Adapter- and quality-trim the FASTQ files using Atria
<a id="footnotes-2"></a>
### Footnotes
<details>
<summary><i>Notes: 2.b. Adapter- and quality-trim the FASTQ files using Atria</i></summary>

<a id="calls-to-unset"></a>
#### Calls to `unset`
- The `unset` command is used to remove or "unset" variables or functions. It's like erasing something from a whiteboard&mdash;once you use `unset`, the variable or function is no longer available in the current session.
- Let's say you have a variable named `my_var` and you want to remove it; to do so, simply type `unset my_var`. After this, `my_var` will no longer hold any value or be recognized by the shell.

<a id="calls-to-typesetdeclare"></a>
#### Calls to `typeset`/`declare`
- `typeset` (also known as `declare`) is used to declare shell variables and give them attributes or set certain properties. It's akin to setting up a box with specific characteristics (like size, color, or label) before putting something into it.
- To declare a new variable with a specific attribute, you use `typeset` followed by options and the variable name. For example, `typeset -i my_num` declares `my_num` as an integer.

<a id="common-options-for-typesetdeclare"></a>
#### Common options for `typeset`/`declare`
1. `-a` (array declaration):
    + Use this to declare a variable as an indexed array.
    + Example: `typeset -a my_array` makes `my_array` an indexed array, where you can store a list of values.
2. `-A` (associative array declaration):
    + This is for declaring an associative array (similar to a dictionary in other languages), where each value is accessed with a unique key.
    + Example: `typeset -A my_dictionary` creates an associative array named `my_dictionary`.
3. `-i` (integer declaration):
    + Declares a variable as an integer. This is useful when you want to ensure that a variable only stores whole numbers.
    + Example: `typeset -i my_num` makes count an integer variable called `my_num`.
4. `-r` (read-only declaration):
    + This makes a variable read-only, meaning once you set its value, it cannot be changed or unset.
    + Example: `typeset -r constant_var=5` creates a read-only variable `constant_var` with the value `5`.

Using `typeset` or `declare` helps you control your variables better. It's like specifying what type of content a folder should have in a filing cabinet, making your scripts more robust and less prone to errors.

There are other options available with `typeset`, such as `-x` for exporting a variable to child processes, or `-f` to list functions. The options available can vary slightly between different shell types (like Bash or Ksh), so it's always a good idea to check the manual (`man bash` or `man ksh`) for the specifics of your environment.

In summary, `unset` is a tool for removing variables or functions, while `typeset`/`declare` is like a Swiss Army knife for creating and managing variables with specific attributes or properties. They are fundamental tools in scripting, enhancing the control and predictability of how scripts behave.

<a id="more-on-operators-particularly-the-logical-operators--and-%7C%7C"></a>
#### More on operators, particularly the logical operators `&&` and `||`
<a id="-operator-logical-and"></a>
##### `&&` Operator (Logical AND)
- The double ampersand `&&` operator allows you to execute a command or set of commands only if the previous command was successful (i.e., it returned an exit status of `0`, which denotes success in \*nix environments).
- `command_1` && `command_2` means "run `command_1`, and if it is successful, then run `command_2`."
- Don't confuse `&&` with `&`: The single ampersand `&` is used to run a command in the background. For example, `command_1 &` runs `command_1` and immediately returns control to the shell, allowing you to continue other work while `command_1` runs in the background.

<a id="%7C%7C-operator-logical-or"></a>
##### `||` Operator (Logical OR)
- The double vertical bar `||` operator lets you execute a command or set of commands only if the previous command failed (i.e., it returned a non-zero exit status).
`command_1 || command_2` means "run `command_1`, and if it fails, then run `command_2`."
- Don't confuse `||` with `|`: The single vertical bar `|` is a pipe, which is used to pass the output of one command as input to another. For example, `command_1 | command_2` takes the output of `command_1` and uses it as input for `command_2`.

<a id="how--and-%7C%7C-work-together"></a>
##### How `&&` and `||` work together
- You can chain these operators for more complex logic. For example, `command_1 && command_2 || command_3` means "run `command_1`, and if it succeeds, then run `command_2`, but if `command_1` fails, then run `command_3`."

<a id="other-relevant-operators--and--logical-not"></a>
##### Other relevant operators: `;` and `!` (logical NOT)
- `;` (semicolon): Used to run commands sequentially, regardless of the success or failure of the previous command. `command_1; command_2` will, no matter what, run `command_1` and then `command_2`.
- `!` (logical NOT): Inverts the exit status of a command. With `!`, if a command fails, then it's associated with `0` (success) exit code, and vice versa: if a command succeeds, then it's associated with a non-zero (failure) exit code.

<a id="best-practices"></a>
##### Best Practices
- Readability: Sometimes, especially for complex logic, it can be clearer to use `if-then-else` statements rather than chaining `&&` and `||`.
- Error handling: Be mindful of how you use these operators in scripts, as they can affect the flow and error handling of your script. Always test how your script behaves in different scenarios.

<a id="summary"></a>
##### Summary
`&&` and `||` are used to control the flow of commands based on their success (`&&`) or failure (`||`). They provide a way to build simple conditional logic directly into the command line. Just remember that they are different from the background operator `&` and the pipe `|`, both of which serve different purposes in scripting.
</details>
<br />
<br />

<a id="3-align-trimmed-fastq-files"></a>
# 3. Align trimmed FASTQ files
<a id="footnotes-3"></a>
<a id="a-generate-a-concatenated-annotated-assembly-of-the-s-cerevisiae-and-s-pombe-genomes"></a>
## a. Generate a concatenated annotated assembly of the *S. cerevisiae* and *S. pombe* genomes
<a id="footnotes-3"></a>
### Footnotes
<details>
<summary><i>Notes: 3.a. Generate a concatenated annotated assembly of the *S. cerevisiae* and *S. pombe* genomes</i></summary>

<a id="breaking-down-the-call-to-mkdir--p-which-makes-use-of-brace-expansion"></a>
#### Breaking down the call to `mkdir -p`, which makes use of brace expansion
1. `[[ ! -d "${dir_genomes}" ]]`: This checks if the directory `"${dir_genomes}"` does not exist. `!` is used for negation, `-d` checks for the existence of a directory, and `${dir_genomes}` is a variable that should hold the path of the directory you're checking.
2. `mkdir -p "${HOME}/genomes/"{"${dir_SP}","${dir_SC}"}/{fasta,gff3}`: This command creates multiple directories in one go.
    + `mkdir -p`: The `mkdir` command is used to create directories, and the `-p` flag ensures that any necessary parent directories are also created (and also prevents an error if the directory already exists).
    + `"${dir_genomes}/"{"${dir_SP}","${dir_SC}"}/{fasta,gff3}`: This is an example of brace expansion. The `mkdir -p` command creates directories in the directory assigned to `dir_genomes`. For each of `${dir_SP}` and `${dir_SC}`, `mkdir -p` will create `fasta` and `gff3` subdirectories. For example, if `${dir_SP}` is "Schizosaccharomyces_pombe" and `${dir_SC}` is "Saccharomyces_cerevisiae", the command will create the following directories:
        - `${HOME}/genomes/Schizosaccharomyces_pombe/fasta`
        - `${HOME}/genomes/Schizosaccharomyces_pombe/gff3`
        - `${HOME}/genomes/Saccharomyces_cerevisiae/fasta`
        - `${HOME}/genomes/Saccharomyces_cerevisiae/gff3`

\#TODO To be continued after modularizing the above
</details>
<br />
<br />
