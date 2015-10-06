---
title: "Permissions and Environment variables"
author: "Kristina Koch, Radhika Khetani"
date: "Tuesday, October 6, 2015"
---

> ## Learning Objectives
> 
> * How to grant or restrict access to files on a multi-user UNIX system
> * What is an "Environment Variable" in a shell.
> * What is $PATH, and why I should care.

Unix controls who can read, modify, and run files using *permissions*.

Let's start with how users are identified in a shared, multi-user system.
We all have a unique username, e.g. rsk27 and a userid 124292.

Find out yours:

```
$ id <username>
```


Users of a multi-user UNIX system can belong to any number of groups,
each of which has a unique group name, and a numeric group ID.

The list of who's in what group is usually stored in the system file `/etc/group`.

Let's see what groups we all belong to:

` $ groups`

Depending on our affiliation, we all belong to at least a couple of groups. I belong to 4 groups,

* rsk27
* bcbio
* hbctraining
* Domain_Users

Unique users and groups are necessary to make sure that Files and Directories that belong to a specific user 

As you can imagine, on a shared system it is important to protect each user's information. To start, every file and directory on a Unix computer belongs to one owner and one group.
Along with each file's content, the operating system stores the numeric IDs of the user and group that own it, which is the "metadata" for a given file.

The user-and-group model means that
for each file
every user on the system falls into one of three categories:
the owner of the file,
someone in the file's group,
and everyone else.

For each of these three categories,
the computer keeps track of
whether people in that category can read the file,
write to the file,
or execute the file
(i.e., run it if it is a program).

For example, if a file had the following set of permissions:

<table class="table table-striped">
<tr><td></td><th>user</th><th>group</th><th>all</th></tr>
<tr><th>read</th><td>yes</td><td>yes</td><td>no</td></tr>
<tr><th>write</th><td>yes</td><td>no</td><td>no</td></tr>
<tr><th>execute</th><td>no</td><td>no</td><td>no</td></tr>
</table>

it would mean that:

*   the file's owner can read and write it, but not run it;
*   other people in the file's group can read it, but not modify it or run it; and
*   everybody else can do nothing with it at all.

Let's look at this model in action.

If we say,

```
$ ls -l /bin/ls
```

It tells us `-rwxr-xr-x. 1 root root 109208 Oct 15  2014 /bin/ls`. 
 
So, `ls` is an executable file that belong to user root and group root, and only they can modify (write) it.


> ## Necessary But Not Sufficient
>
> The fact that something is marked as executable
> doesn't actually mean it contains or is a program of some kind.
> We could easily mark the `~/unix_oct2015/raw_fastq/Irrel_kd_1.subset.fq` file as executable
> using the commands that are introduced below.
> Depending on the operating system we're using,
> trying to "run" it will fail
> (because it doesn't contain instructions the computer recognizes).


Now let's run the command `ls -l ~/unix_oct2015`, to list the files in that directory:

```
$ ls -l

drwxrwsr-x 2 rsk27 rsk27  78 Oct  6 10:29 genomics_data
drwxrwsr-x 2 rsk27 rsk27 228 Oct  6 10:28 raw_fastq
-rw-rw-r-- 1 rsk27 rsk27 377 Oct  6 10:28 README.txt
drwxrwsr-x 2 rsk27 rsk27 238 Oct  6 10:28 reference_data
drwxrwsr-x 2 rsk27 rsk27 555 Oct  6 10:29 reference_STAR
```

The `-l` flag tells `ls` to give us a long-form listing.
It's a lot of information, so let's go through the columns in turn.

||On the right side, we have the files' names. Next to them,
moving left, are the times and dates they were last modified. Backup systems and other tools use this information in a variety of ways, but you can use it to tell when you (or anyone else with permission) last changed a file.

Next to the modification time is the file's size in bytes
and the names of the user and group that owns it
(in this case, `rsk27` and `rsk27` respectively).
We'll skip over the second column for now
(the one showing `1` for each file), 
because it's the first column that we care about most.
This shows the file's permissions, i.e., who can read, write, or execute it.

Let's have a closer look at one of those permission strings for README.txt:
`-rw-rw-r--`.
The first character tells us what type of thing this is:
'-' means it's a regular file,
while 'd' means it's a directory,
and other characters mean more esoteric things.

The next three characters tell us what permissions the file's owner has.
Here, the owner can read and write the file: `rw-`.

The middle triplet shows us the group's permissions.
If the permission is turned off, we see a dash, so `rw-` means "read and write, but not execute". (In this case the group and the owner are the same so it makes sense that this is the same for both categories.)

The final triplet shows us what everyone who isn't the file's owner, or in the file's group, can do. In this case, it's `r--` again, so everyone on the system can look at the file's contents.

To change permissions, we use the `chmod` command (whose name stands for "change mode").
Let's make our README.txt file inaccessible to all users on the system, they can currently read it:

```
$ ls -l ~/unix_oct2015/README.txt

-rw-rw-r-- 1 rsk27 rsk27 377 Oct  6 10:28 /home/rsk27/unix_oct2015/README.txt

$ chmod o-rw ~/unix_oct2015/README.txt

$ ls -l ~/unix_oct2015/README.txt

-rw-rw---- 1 rsk27 rsk27 377 Oct  6 10:28 /home/rsk27/unix_oct2015/README.txt
```

The 'o' signals that we're changing the privileges of "others".

Let's change it back to allow it to be readable by others:

```
$ chmod o+r ~/unix_oct2015/README.txt

$ ls -l ~/unix_oct2015/README.txt

-rw-rw-r-- 1 rsk27 rsk27 377 Oct  6 10:28 /home/rsk27/unix_oct2015/README.txt
```

If we wanted to make this an executable file for ourselves, the file's owners, we would say `chmod u+rwx`, where the 'u' signals that we are changing permission for the file's owner. To change permissions for a whole group, you'd use the letter "g" `chmod g-w`. 

Before we go any further,
let's run `ls -a -l` on the `~/unix_oct2015` directory to get a long-form listing:

```
$ ls -l

drwxrwsr-x 2 rsk27 rsk27  78 Oct  6 10:29 genomics_data
drwxrwsr-x 2 rsk27 rsk27 228 Oct  6 10:28 raw_fastq
-rw-rw-r-- 1 rsk27 rsk27 377 Oct  6 10:28 README.txt
drwxrwsr-x 2 rsk27 rsk27 238 Oct  6 10:28 reference_data
drwxrwsr-x 2 rsk27 rsk27 555 Oct  6 10:29 reference_STAR
```

Look at the permissions for directories (`drwxrwsr-x`): the 'x' indicates that "execute" is turned on. What does that mean? A directory isn't a program or an executable file, we can't "run" it.

Well, 'x' means something different for directories.
It gives someone the right to *traverse* the directory, but not to look at its contents.
The distinction is subtle, so let's have a look at an example.

Dr. Vlad Smith's home directory has three subdirectories called `venus`, `mars`, and `pluto`:

![execute](../img/permission-directory.png "Execute Permission for Directories")

Each of these has a subdirectory in turn called `notes`, and those sub-subdirectories contain various files.
If a user's permissions on `venus` are 'r-x', then if she tries to see the contents of `venus` and `venus/notes` using `ls`, the computer lets her see both.
If her permissions on `mars` are just 'r--', then she is allowed to read the contents of both `mars` and `mars/notes`.
But if her permissions on `pluto` are only '--x', she cannot see what's in the `pluto` directory: `ls pluto` will tell her she doesn't have permission to view its contents.
If she tries to look in `pluto/notes`, though, the computer will let her do that.
She's allowed to go through `pluto`, but not to look at what's there. She will be able to do this, only if she knows that there is a file called `notes` in the directory, since she cannot list what is in there.

This trick gives people a way to make some of their directories visible to the world as a whole without opening up everything else.


> ## Challenge
> If `ls -l myfile.php` returns the following details:
>
> ```
> -rwxr-xr-- 1 caro zoo  2312  2014-10-25 18:30 myfile.php
> ```
> 
> Which of the following statements is true?
> 
> 1. caro (the owner) can read, write, and execute myfile.php
> 2. caro (the owner) cannot write to myfile.php
> 3. members of caro (a group) can read, write, and execute myfile.php
> 4. members of zoo (a group) cannot execute myfile.php

