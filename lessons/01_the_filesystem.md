---
layout: page
title: "The Shell"
comments: true
date: 2015-10-4
---

# The Shell
Author: Sheldon  McKay, Mary Piper 

Adapted from the lesson by Tracy Teal.
Original contributors:
Paul Wilson, Milad Fatenejad, Sasha Wood and Radhika Khetani for Software Carpentry (http://software-carpentry.org/)

## Learning Objectives
- What is the shell?
- How do you access it?
- How do you use it?
  - Getting around the Unix file system
  - looking at files
  - manipulating files
  - automating tasks
- What is it good for?
- Where are resources where I can learn more? (because the shell is awesome)

## What is the shell?

The *shell* is a program that presents a command line interface which allows you to control your computer using commands entered with a keyboard instead of controlling graphical user interfaces (GUIs) with a mouse/keyboard combination.

There are many reasons to learn about the shell.

* For most bioinformatics tools, you have to use the shell. There is no
graphical interface. If you want to work in metagenomics or genomics you're
going to need to use the shell.
* The shell gives you *power*. The command line gives you the power to do your work more efficiently and
more quickly.  When you need to do things tens to hundreds of times,
knowing how to use the shell is transformative.
* To use remote computers or cloud computing, you need to use the shell.


![Automation](../img/gvng.jpg)

   Unix is user-friendly. It's just very selective about who its friends are.


Today we're going to go through how to access Unix/Linux and some of the basic
shell commands.

## Information on the shell

shell cheat sheets:<br>
* [http://fosswire.com/post/2007/08/unixlinux-command-cheat-sheet/](http://fosswire.com/post/2007/08/unixlinux-command-cheat-sheet/)
* [https://github.com/swcarpentry/boot-camps/blob/master/shell/shell_cheatsheet.md](https://github.com/swcarpentry/boot-camps/blob/master/shell/shell_cheatsheet.md)

Explain shell - a web site where you can see what the different components of
a shell command are doing.  
* [http://explainshell.com](http://explainshell.com)
* [http://www.commandlinefu.com](http://www.commandlinefu.com)


## Setting up

We will spend most of our time learning about the basics of the shell by manipulating some experimental data.

Since we are going to be working with this data on our remote server, **Orchestra**, we first need to log onto the server. After we're logged on, we will each make our own copy of the example data folder.

### Logging onto Orchestra with Macs

Using the Terminal, you can use the command 'ssh' and your eCommons username to login. Type:

```ssh username@orchestra.med.harvard.edu```

You will receive a prompt for your password, and you should type in your ECommons password. 

### Logging onto Orchestra with Windows

By default, there is no terminal for the bash shell available in the Windows OS, so you have to use a downloaded program, **[Putty](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html)**.

When you open Putty, you will see the following GUI.

![Putty window](../img/putty-1.PNG)

Type in "orchestra.med.harvard.edu" in the window under "Host Name (or IP address) and click on "Open"

![Connect to Orchestra](../img/putty-2.PNG)

A warning window will pop up the first time you try to connect to a cluster (remote server), say "Yes". Once you do that, you should be able to enter your login ID which is your eCommons ID. Add ID and press enter.

![Log in](../img/putty-5.PNG)

Once you press enter, it will prompt you for a password. Type in your password, when you do this nothing will appear on the screen until you press enter. When you press enter, the interface will change and you have started a bash terminal.


### Copying example data folder

Once logged in, you should see the Orchestra news and the command prompt: 

```$ ```

The command prompt will have some characters before it, something like "-bash-4.1", this is telling you what the name of the computer you are working on is.

The first command we will type on the command prompt will be to start a so-called "interactive session" on Orchestra.

```$ bsub -Is -q interactive bash```

Press enter after you type in that command. You will get a couple of messages, but in a few seconds you should get back the command prompt `$`; the string of characters before the command prompt, however, have changed. They should say something like `rsk27@clarinet002-062`. *We will be explaining what this means in more detail later this afternoon when we talk about HPC and Orchestra.* 

Make sure that your command prompt is now preceded by a character string that contain words like "clarinet", "bassoon", etc.

Copy our example data folder to your home directory using the following command:

```$ cp -r /groups/hbctraining/unix_oct2015/ .```

>'cp' is the command for copy. This command required you to specify the location of the item you want to copy (/groups/hbctraining/unix_oct2015/) and the location of the destination (.) please note the space between the 2 in the command. The "-r" is an option that modifies the copy command to do something slightly different than usual. The "." means "here", i.e. the destination location is where you currently are.

## Starting with the shell

We have each created our own copy of the example data folder into our home directory, **unix_oct2015**. Let's go into the data folder and explore the data using the shell.

```$ cd ~/unix_oct2015```

'cd' stands for 'change directory'

Let's see what is in here. Type:

```$ ls```

You will see:

	genomics_data  other  raw_fastq  README.txt  reference_data

ls stands for 'list' and it lists the contents of a directory.

There are five items listed.  What are they? We can use a command line argument with `ls` to get more information.

```$ ls -F```

	genomics_data/  other/  raw_fastq/  README.txt  reference_data/

Anything with a "/" after it is a directory. Things with a "*" after them are programs.  If there are no decorations, it's a file.

You can also use the command

```$ ls -l```

	total 124
	drwxrwsr-x 2 mp298 mp298  78 Sep 30 10:47 genomics_data
	drwxrwsr-x 6 mp298 mp298 107 Sep 30 10:47 other
	drwxrwsr-x 2 mp298 mp298 228 Sep 30 10:47 raw_fastq
	-rw-rw-r-- 1 mp298 mp298 377 Sep 30 10:47 README.txt
	drwxrwsr-x 2 mp298 mp298 238 Sep 30 10:47 reference_data

to see whether items in a directory are files or directories. `ls -l` gives a lot more
information too.

Let's go into the raw_fastq directory and see what is in there.

```$ cd raw_fastq/```

```$ ls -F```

	Irrel_kd_1.subset.fq  Irrel_kd_3.subset.fq  Mov10_oe_2.subset.fq
	Irrel_kd_2.subset.fq  Mov10_oe_1.subset.fq  Mov10_oe_3.subset.fq

All six items in this directory have no trailing slashes, so they are all files.


## Arguments

Most programs take additional arguments that control their exact behavior. For example, `-F` and `-l` are arguments to `ls`.  The `ls` program, like many programs, take a lot of arguments. Another useful one is '-a', which show everything, including hidden files.  How do we know what the options are to particular commands?

Most commonly used shell programs have a manual. You can access the
manual using the `man` program. Try entering:

```$ man ls```

This will open the manual page for `ls`. Use the space key to go forward and b to go backwards. When you are done reading, just hit `q` to quit.

Programs that are run from the shell can get extremely complicated. To see an example, open up the manual page for the `find` program. No one can possibly learn all of these arguments, of course. So you will probably find yourself referring back to the manual page frequently.


## The Unix directory file structure (a.k.a. where am I?)

As you've already just seen, you can move around in different directories or folders at the command line. Why would you want to do this, rather than just navigating around the normal way.

When you're working with bioinformatics programs, you're working with your data and it's key to be able to have that data in the right place and make sure the program has access to the data. Many of the problems people run in to with command line bioinformatics programs is not having the data in the place the program expects it to be.


## Moving around the file system

Let's practice moving around a bit.

We're going to work in that `unix_oct2015` directory.

First we did something like go to the folder of our username. Then we opened
`unix_oct2015` then `raw_fastq`

Let's draw out how that went.

Now let's draw some of the other files and folders we could have clicked on.

This is called a hierarchical file system structure, like an upside down tree
with root (/) at the base that looks like this.

![Unix](../img/Slide1.jpg)

That (/) at the base is often also called the 'top' level.

When you are working at your computer or log in to a remote computer,
you are on one of the branches of that tree, your home directory (e.g. /home/username)

Now let's go do that same navigation at the command line.

Type

```$ cd```

This puts you in your home directory. This folder here.

Now using `cd` and `ls`, go in to the `unix_oct2015` directory and list its contents. Now go into the `raw_fastq` and list its contents.

Let's also check to see where we are. Sometimes when we're wandering around in the file system, it's easy to lose track of where we are and get lost.

If you want to know what directory you're currently in, type

```$ pwd```

This stands for 'print working directory'. The directory you're currently working in.

What if we want to move back up and out of the `raw_fastq` directory? Can we just type `cd unix_oct2015`? Try it and see what happens.

To go 'back up a level' we need to use `..`

Type

```$ cd ..```

Now do `ls` and `pwd`. See now that we went back up in to the `unix_oct2015`
directory. `..` means go back up a level.

* * * *
**Exercise**

Now we're going to try a hunt. Find a hidden directory in `unix_oct2015` list its contents, and find the text file in there.  What is the name of the file?

Hint: hidden files and folders in unix start with `.`, for example `.my_hidden_directory`
* * * *


## Examining the contents of other directories

By default, the `ls` commands lists the contents of the working directory (i.e. the directory you are in). You can always find the directory you are in using the `pwd` command. However, you can also give `ls` the names of other directories to view. Navigate to the home directory if you are not already there.

Type:

```$ cd```

Then enter the command:

```$ ls unix_oct2015/```

This will list the contents of the `unix_oct2015` directory without you having to navigate there.

The `cd` command works in a similar way. Try entering:

```$ cd```

```$ cd unix_oct2015/raw_fastq/```

and you will jump directly to `raw_fastq` without having to go through the intermediate directory.

****
**Exercise**

List the `Mov10_oe_1.subset.fq` file from your home directory without changing directories
****

### Shortcut: Tab Completion

Navigate to the home directory. Typing out directory names can waste a lot of time. When you start typing out the name of a directory, then hit the tab key, the shell will try to fill in the rest of the directory name. For example, type `cd` to get back to your home directly, then enter:

```$ cd uni<tab>```

The shell will fill in the rest of the directory name for `unix_oct2015`. Now go to `unix_oct2015/raw_fastq` and 

```$ ls Mov10_oe_<tab><tab>```

When you hit the first tab, nothing happens. The reason is that there are multiple directories in the home directory which start with `Mov10_oe_`. Thus, the shell does not know which one to fill in. When you hit tab again, the shell will list the possible choices.

Tab completion can also fill in the names of programs. For example, enter `e<tab><tab>`. You will see the name of every program that starts with an `e`. One of those is `echo`. If you enter `ec<tab>` you will see that tab completion works.

## Full vs. Relative Paths

The `cd` command takes an argument which is the directory name. Directories can be specified using either a *relative* path or a full *path*. The directories on the computer are arranged into a hierarchy. The full path tells you where a directory is in that hierarchy. Navigate to the home directory (`cd`). Now, enter the `pwd` command and you should see:

`$ pwd`

`/home/username`

which is the full name of your home directory. This tells you that you are in a directory called `username`, which sits inside a directory called `home` which sits inside the very top directory in the hierarchy. The very top of the hierarchy is a directory called `/` which is usually referred to as the *root directory*. So, to summarize: `username` is a directory in `home` which is a directory in `/`.

Now enter the following command:

```$ cd /home/username/unix_oct2015/.my_hidden_directory```

This jumps to `.my_hidden_directory`. Now go back to the home directory (`cd`). We saw
earlier that the command:

```$ cd unix_oct2015/.my_hidden_directory/```

had the same effect - it took us to the `my_hidden_directory` directory. But, instead of specifying the full path (`/home/username/unix_oct2015/.my_hidden_directory`), we specified a *relative path*. In other words, we specified the path relative to our current directory. A full path always starts with a `/`. A relative path does not.

A relative path is like getting directions from someone on the street. They tell you to "go right at the Stop sign, and then turn left on Main Street". That works great if you're standing there together, but not so well if you're trying to tell someone how to get there from another country. A full path is like GPS coordinates. It tells you exactly where something is no matter where you are right now.

You can usually use either a full path or a relative path depending on what is most convenient. If we are in the home directory, it is more convenient to just enter the relative path since it involves less typing.

Over time, it will become easier for you to keep a mental note of the structure of the directories that you are using and how to quickly navigate amongst them.

***
**Exercise**

Now, list the contents of the `/bin` directory. Do you see anything familiar in there? 
How can you tell these are programs rather than plain files?
***

## Saving time with shortcuts, wild cards, and tab completion

### Shortcuts

There are some shortcuts which you should know about. Dealing with the
home directory is very common. So, in the shell the tilde character,
"~", is a shortcut for your home directory. Navigate to the `raw_fastq`
directory:

```$ cd```

```$ cd unix_oct2015/raw_fastq```

Then enter the command:

```$ ls ~```

This prints the contents of your home directory, without you having to type the full path. The shortcut `..` always refers to the directory above your current directory. Thus:

```$ ls ..```

prints the contents of the `/home/username/unix_oct2015`. You can chain these together, so:

```$ ls ../..```

prints the contents of `/home/username` which is your home directory. Finally, the special directory `.` always refers to your current directory. So, `ls`, `ls .`, and `ls ././././.` all do the same thing, they print the contents of the current directory. This may seem like a useless shortcut right now, but we'll see when it is needed in a little while.

To summarize, while you are in your home directory, the commands `ls ~`, `ls ~/.`, `ls ../../`, and `ls /home/username` all do exactly the same thing. These shortcuts are not necessary, they are provided for your convenience.

### Our data set: FASTQ files

We did an experiment and want to look at sequencing results. We want to be able to look at these files and do some things with them.


### Wild cards

Navigate to the `~/unix_oct2015/raw_fastq` directory. This
directory contains our FASTQ files.

The '*' character is a shortcut for "everything". Thus, if you enter `ls *`, you will see all of the contents of a given directory. Now try this command:

```$ ls *fq```

This lists every file that ends with a `fq`. This command:

```$ ls /usr/bin/*.sh```

Lists every file in `/usr/bin` that ends in the characters `.sh`.

```$ ls Mov10*fq```

lists only the files that begin with 'Mov10' and end with 'fq'

So how does this actually work? Well...when the shell (bash) sees a word that contains the `*` character, it automatically looks for filenames that match the given pattern. 

We can use the command `echo` to see wilcards are they are intepreted by the shell.

```$ echo *.fq```

	Irrel_kd_1.subset.fq Irrel_kd_2.subset.fq Irrel_kd_3.subset.fq
	Mov10_oe_1.subset.fq Mov10_oe_2.subset.fq Mov10_oe_3.subset.fq

The '*' is expanded to include any file that ends with '.fq'


****
**Exercise**

Do each of the following using a single `ls` command without
navigating to a different directory.

1.  List all of the files in `/bin` that start with the letter 'c'
2.  List all of the files in `/bin` that contain the letter 'a'
3.  List all of the files in `/bin` that end with the letter 'o'

BONUS: List all of the files in `/bin` that contain the letter 'a' or 'c'

****


## Command History

You can easily access previous commands.  Hit the up arrow. Hit it again.  You can step backwards through your command history. The down arrow takes your forwards in the command history.

\^-C will cancel the command you are writing, and give you a fresh prompt.

\^-R will do a reverse-search through your command history.  This
is very useful.

You can also review your recent commands with the `history` command.  Just enter:

```$ history```

to see a numbered list of recent commands, including this just issues
`history` command.  You can reuse one of these commands directly by
referring to the number of that command.

If your history looked like this:

    259  ls *
    260  ls /usr/bin/*.sh
    261  ls *fq

then you could repeat command #260 by simply entering:

```$ !260```


(that's an exclamation mark).  You will be glad you learned this when you try to re-run very complicated commands.

****
**Exercise**

1. Find the line number in your history for the last exercise (listing
files in `/bin`) and reissue that command.

****



## Examining Files

We now know how to switch directories, run programs, and look at the
contents of directories, but how do we look at the contents of files?

The easiest way to examine a file is to just print out all of the
contents using the program `cat`. Change directories to the `unix_oct2015/other/` and enter the following command:

`$ cat sequences.fa`

This prints out the all the contents of `sequences.fa` to the screen.

What does this file contain?

`cat` is a terrific program, but when the file is really big, it can be annoying to use. The program, `less`, is useful for this case. Let's take a look at the raw_fastq files. These files are quite large, so we probably do not want to use the `cat` command to look at them. Instead, we can use the `less` command. 

Move back to our `raw_fastq` directory and enter the following command:

`less Mov10_oe_1.subset.fq`

We will explore fastq files in more detail later, but notice that fastq files have four lines of data associated with every sequence read. Not only is there a header line and the nucleotide sequence, similar to a fasta file, but fastq files also contain quality information for each nucleotide in the sequence. 

The `less` command opens the file, and lets you navigate through it. The commands are identical to the `man` program.

**Some commands in `less`**

| key     | action |
| ------- | ---------- |
| "space" | to go forward |
|  "b"    | to go backwarsd |
|  "g"    | to go to the beginning |
|  "G"    | to go to the end |
|  "q"    | to quit |

`less` also gives you a way of searching through files. Just hit the "/" key to begin a search. Enter the name of the word you would like to search for and hit enter. It will jump to the next location where that word is found. If you hit "/" then "enter", `less` will just repeat the previous search. `less` searches from the current location and works its way forward. If you are at the end of the file and search for the word "cat", `less` will not find it. You need to go to the beginning of the file and search.

For instance, let's search for the sequence `GAGACCCCACGGGAGGCCA` in our file. You can see that we go right to that sequence and can see what it looks like. (Remember to hit 'q' to exit the `less` program)

Remember, the `man` program actually uses `less` internally and
therefore uses the same commands, so you can search documentation
using "/" as well!

There's another way that we can look at files, and in this case, just
look at part of them. This can be particularly useful if we just want
to see the beginning or end of the file, or see how it's formatted.

The commands are `head` and `tail` and they just let you look at
the beginning and end of a file respectively.

```$ head Mov10_oe_1.subset.fq ```

```$ tail Mov10_oe_1.subset.fq```

The `-n` option to either of these commands can be used to print the first or last `n` lines of a file. To print the first/last line of the file use:

```$ head -n 1 Mov10_oe_1.subset.fq```

```$ tail -n 1 Mov10_oe_1.subset.fq```

## Creating, moving, copying, and removing

Now we can move around in the file structure, look at files, search files, redirect. But what if we want to do normal things like copy files or move them around or get rid of them. Sure we could do most of these things without the command line, but what fun would that be?! Besides it's often faster to do it at the command line, or you'll be on a remote server like Amazon where you won't have another option.


Our raw data in this case is fastq files.  We don't want to change the original files,
so let's make a copy to work with.

Lets copy the file using the `cp` command. The `cp` command backs up the file. Navigate to the `data` directory and enter:

```$ cp Mov10_oe_1.subset.fq Mov10_oe_1.subset-copy.fq```

```$ ls -F```

Now ``Mov10_oe_1.subset-copy.fq`` has been created as a copy of `Mov10_oe_1.subset.fq`

Let's make a 'backup' directory where we can put this file.

The `mkdir` command is used to make a directory. Just enter `mkdir`
followed by a space, then the directory name.

```$ mkdir backup```

We can now move our backed up file in to this directory. We can move files around using the command `mv`. Enter this command:

```$ mv *-copy.fq backup```

```$ ls -al backup```

	drwxrwsr-x 2 mp298 mp298       43 Sep 30 13:59 .
	drwxrwsr-x 8 mp298 mp298      203 Sep 30 13:58 ..
	-rw-rw-r-- 1 mp298 mp298 75706556 Sep 30 13:56 Mov10_oe_1.subset-copy.fq

The `mv` command is also how you rename files. Since this file is so
important, let's rename it:

```$ cd backup```

`$ mv Mov10_oe_1.subset-copy.fq Mov10_oe_1.subset-copy.fq_DO_NOT_TOUCH!`

`$ ls`

	Mov10_oe_1.subset-copy.fq_DO_NOT_TOUCH!

Finally, we decided this was silly and want to start over.

```$ cd ..```

```$ rm backup/Mov*```

The `rm` file permanently removes the file. Be careful with this command. It doesn't
just nicely put the files in the Trash. They're really gone.



* * * *
**Exercise**

Do the following:

1.  Create a backup of your fastq files
2.  Create a backup directory 
3.  Copy your backup files there

* * * *

By default, `rm`, will NOT delete directories. You can tell `rm` to delete a directory using the `-r` option. Let's delete that `new` directory we just made. Enter the following command:

```$ rm -r backup/```

## Writing files

We've been able to do a lot of work with files that already exist, but what if we want to write our own files. Obviously, we're not going to type in a FASTA file, but you'll see as we go through other tutorials, there are a lot of reasons we'll want to write a file, or edit an existing file.

To write in files, we're going to use the program `nano`. We're going to create
a file that contains the favorite grep command so you can remember it for later. We'll name this file 'awesome.sh'.

```$ nano awesome.sh```

Now you have something that looks like

![nano1.png](../img/nano1.png)

Type in your command:

```bash grep -A 3 -B 1 GAGACCCCACGGGAGGCCA Mov10_oe_1.subset.fq```

 so it looks like

![nano2.png](../img/nano2.png)

Now we want to save the file and exit. At the bottom of nano, you see the "\^X Exit". That means that we use Ctrl-X to exit. Type `Ctrl-X`. It will ask if you want to save it. Type `y` for yes. Then it asks if you want that file name. Hit 'Enter'.

Now you've written a file. You can take a look at it with less or cat, or open it up again and edit it.

***
**Exercise**

Open 'awesome.sh' and add "echo AWESOME!" after the grep command and save the file.

We're going to come back and use this file in just a bit.

***

**Commands, options, and keystrokes covered in this lesson**

```bash
cd
ls
man
pwd
~ (home dir)
. (current dir)
.. (parent dir)
*  (wildcard)
echo
ctrl-C (cancel current command)
ctrl-R (reverse history search)
ctrl-A (start of line)
ctrl-E (end of line)
history
! (repeat cmd)
cat
less
head
tail
cp
mdkir
mv
rm
nano
```


