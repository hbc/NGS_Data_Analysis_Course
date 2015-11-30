---
title: "The Shell: Loops & Scripts"
author: "Bob Freeman, Mary Piper, Radhika Khetani"
date: "Wednesday, October 28, 2015"
---

Approximate time: 60 minutes

## Learning Objectives

* Learn how to operate on multiple files 
* Capture previous commands into a script to re-run later
* Abstract your script for flexibility
* Write a series of scripts, that are increasingly more flexible, to automate your workflow

Now that you've been using quite a number of commands to interrogate your data, 
wouldn't it be great if you could do this for each set of data that comes in, without having to manually re-type the commands?

Welcome to the beauty and purpose of shell scripts.

## Shell scripts

We are finally ready to see what makes the shell such a powerful programming environment. We are going to take the commands we repeat frequently and save them in files so that we can re-run all those operations again later by typing one single command. For historical reasons, a bunch of commands saved in a file is usually called a shell script, but make no mistake: this is actually a small program.

Shell scripts are text files that contain commands we want to run over and over again. As with any file, you can give a shell script any name and usually have the extension `.sh`. Let's write a shell script that tells us what our current working directory is and then lists the contents of the directory. First open a new file using nano:

	$ nano listing.sh
	
Then type in the following lines in the `listing.sh` file:

	echo "Your current working directory is:"
	pwd
	echo "These are the contents of this directory:"
	ls -l 

Close nano and save the file. Now let's run the new script we have created. To run a shell script you usually use the `bash` or `sh` command.

	$ sh listing.sh
	
> Did it work like you expected?
> 
> Were the `echo` commands helpful in letting you know what came next?

This is a very simple shell script. In this session and in upcoming sessions, we will be learning how to write more complex ones. You will see how the power of scripts can make our lives much easier.

## Bash variables
A *variable* is a common concept shared by many programming languages. Variables are essentially a symbolic/temporary name for, or reference to, information. Variables are analogous to "buckets", where information can be stored, maintained and modified without too much hassle. 

Extending the bucket analogy: the bucket has a name associated with it, i.e. the name of the variable, and when referring to the information in the bucket, we use the name of the bucket, and do not directly refer to the actual data stored in it.

In the example below, we define a variable or a 'bucket' called `file`. We will put a filename `Mov10_oe_1.subset.fq` as the value inside the bucket.

	$ file=Mov10_oe_1.subset.fq

Once you press return, you should be back at the command prompt. *How do we know that we actually created the bash variable?* We can use the echo command to list what's inside `file`:

	$ echo $file

What do you see in the terminal? If the variable was not created, the command will return nothing. Did you notice that when we created the variable we just typed in the variable name, but when using it as an argument to the `echo` command, we explicitly use a `$` in front of it (`$file`). Why? 

Well, in the former, we're setting the value, while in the latter, we're retrieving the value. This is standard shell notation (syntax) for defining and using variables. **Don't forget the `$` when you want to retrieve the value of a variable!** 

Let's try another command using the variable that we have created. In the last lesson, we introduced the `wc -l` command which allows us to count the number of lines in a file. We can count the number of lines in `Mov10_oe_1.subset.fq` by referencing the `file` variable, but first move into the `raw_fastq` directory:

	$ cd ~/unix_oct2015/raw_fastq
	$ wc -l $file

Ok, so we know variables are like buckets, and so far we have seen that bucket filled with a single value. **Variables can store more than just a single value.** They can store multiple values and in this way can be useful to carry out many things at once. Let's create a new variable called `filenames` and this time we will store *all of the filenames* in the `raw_fastq` directory as values. 

To list all the filenames in the directory that have a `.fq` extension, we know the command is:

	$ ls *.fq
	
Now we want to re-direct the output of `ls` into a variable. We will give that variable the name `filenames`:

	$ filenames=`ls *.fq`

Check and see what's stored inside our newly created variable using `echo`:
	
	$ echo $filenames

Let's try the `wc -l` command again, but this time using our new variable `filenames` as the argument:

	$ wc -l $filenames
	
What just happened? Because our variable contains multiple values, the shell runs the command on each value stored in `filenames` and prints the results to screen. 

> Try using some of the other commands we learned in previous lessons (i.e `head`, `tail`) on the `filename` variable. 

## Loops

Another powerful concept in the Unix shell is the concept of "Loops". We have just shown you that you can run a single command on multiple files by creating a variable whose values are the filenames that you wish to work on. But what if you want to **run a sequence of multiple commands, on multiple files**? This is where loop come in handy!

Looping is a concept shared by several programming languages, and its implementation in bash is very similar to other languages. 

The structure or the syntax of (*for*) loops in bash is as follows:

```
$ for (variable_name) in (list)
> do
>   (command $variable_name) 
> done
```

where the ***variable_name*** defines (or initializes) a variable that takes the value of every member of the specified ***list*** one at a time. At each iteration, the loop retrieves the value stored in the variable (which is a member of the input list) and runs through the commands indicated between the `do` and `done` one at a time. *This syntax/structure is virtually set in stone.* 

For example, we can run the same commands (`echo` and `wc -l`) used in the "Bash variables" section but this time run them sequentially on each file:

```
$ ls  *.fq		# list all files ending in .fq

$ for var in *.fq
> do
>   echo $var
>   wc -l $var
> done
```

####What does this loop do? 
Most simply, it writes to the terminal (`echo`) the name of the file and the number of lines (`wc -l`) for each files that end in `.fq` in the current directory. The output is almost identical to what we had before.

In this case the list of files is specified using the asterisk wildcard: `*.fq`, i.e. all files that end in `.fq`. Then, we execute 2 commands between the `do` and `done`. With a loop, we execute these commands for each file at a time. Once the commands are executed for one file, the loop then executes the same commands on the next file in the list. 

Essentially, **the number of loops == the number of items in the list**, in our case that is 6 times since we have 6 files in `~/unix_oct2015/raw_fastq` that end in `.fq`. This is done by changing the value of the `var` variable 6 times. 

Of course, `var` is a useless variable name. But it doesn't matter what variable name we use and we can make it something more intuitive.

```bash
$ for filename in *.fq
> do
>   echo $filename
>   wc -l $filename
> done
```
In the long run, it's best to use a name that will help point out its function, so your future self will understand what you are thinking now.

Now that we understand the concept of looping, let's put that to work:

```bash
$ for filename in *.fq
> do 
>   echo $filename >> ../other/number-of-badreads.txt
>   grep -B1 -A2 NNNNNNNNNN $filename | wc -l >> ../other/number-of-badreads.txt
>   grep -B1 -A2 NNNNNNNNNN $filename > ../other/$filename.badreads.fastq 
> done
```
We have used the for loop along with the `>>` redirection symbol to populate one file with all the bad reads in our data set.

Pretty simple and cool, huh?

## Automating with Scripts

Now that you've learned how to use loops and variables, let's put this processing power to work. Imagine, if you will, a series of commands that would do the following for us each time we get a new data set:

- Use for loop to iterate over each FASTQ file
- Dump out bad reads into a new file
- Get the count of the number of bad reads
- And after all the FASTQ files are processed, we generate one summary file of the bad read counts

You might not realize this, but this is something that you now know how to do. Let's get started...

Move to our sample data directory and use `nano` to create our new script file:

`$ cd ~/unix_oct2015/raw_fastq`

`$ nano generate_bad_reads_summary.sh`

We always want to start our scripts with a shebang line: 

`#!/bin/bash`

This line is the absolute path to the Bash interpreter. The shebang line ensures that the bash shell interprets the script even if it is executed using a different shell.

After the shebang line, we enter the commands we want to execute. First we want to move into our `raw_fastq` directory:

```
$ cd ~/unix_oct2015/raw_fastq
```

And now we loop over all the FASTQs:

```bash
for filename in ~/unix_oct2015/raw_fastq/*.fq;
```

and we execute the commands for each loop:

```bash
do
  # tell us what file we're working on
  echo $filename;
  
  # grab all the bad read records into new file
  grep -B1 -A2 NNNNNNNNNN $filename > $filename-badreads.fastq;
``` 
  
We'll also get the counts of these reads and put that in a new file, using the count flag of `grep`:

```bash
# grab the # of bad reads from our badreads file
grep -cH NNNNNNNNNN $filename-badreads.fastq > $filename-badreads.counts;
done
```

If you've noticed, we slipped a new `grep` flag `-H` in there. This flag will report the filename along with the match string. This won't matter within each file, but it will matter when we generate the summary:

```bash
# grab all our bad read info to a summary file
cat *.counts > bad-reads.count.summary
```

And now, as a best practice of capturing all of our work into a running summary log:

```bash
# and add this summary to our run log
cat bad-reads.count.summary >> ../runlog.txt
```

You're script should look like:

```bash
#!/bin/bash

cd ~/unix_oct2015/raw_fastq

for filename in ~/unix_oct2015/raw_fastq/*.fq; do 
echo $filename;
grep -B1 -A2 NNNNNNNNNN $filename > $filename-badreads.fastq;
grep -cH NNNNNNNNNN $filename-badreads.fastq > $filename-badreads.counts;
done

cat *.counts > bad-reads.count.summary

cat bad-reads.count.summary >> ../runlog.txt

```

Exit out of `nano`, and voila! You now have a script you can use to assess the quality of all your new datasets. Your finished script, complete with comments, should look like the following:

```bash
#!/bin/bash 

# enter directory with raw FASTQs
cd ~/unix_oct2015/raw_fastq

# count bad reads for each FASTQ file in our directory
for filename in ~/unix_oct2015/raw_fastq/*.fq; do 
  echo $filename; 

  # grab all the bad read records
  grep -B1 -A2 NNNNNNNNNN $filename > $filename-badreads.fastq;

  # grab the # of bad reads from our bad reads file
  grep -cH NNNNNNNNNN $filename-badreads.fastq > $filename-badreads.counts;
done

# add all our bad read info to a summary file
cat *.counts > bad-reads.count.summary

# and add this summary to our run log
cat bad-reads.count.summary >> ../runlog.tx

```

To run this script, we simply enter the following command:

```bash
$ bash generate_bad_reads_summary.sh
```

To keep your data organized, let's move all of the bad read files out of our `raw_fastq` directory into the `other` directory

`$ mv ~/unix_oct2015/raw_fastq/*bad* ~/unix_oct2015/other`


---
*The materials used in this lesson was derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*

