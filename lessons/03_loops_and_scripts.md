---
title: "The Shell: Loops & Scripts"
author: "Bob Freeman, Mary Piper"
date: "Wednesday, October 7, 2015"
---

Approximate time: 60 minutes

## Learning Objectives

* Learn how to operate on multiple files 
* Capture previous commands into a script to re-run later
* Abstract your script for flexibility
* Write a series of scripts, that are increasingly more flexible, to automate your workflow

Now that you've been using quite a number of commands to interrogate your data, 
wouldn't it be great if you could generate a report for your data and the special
SRARunTable files as a part of the interrogation process? Even better, could we do this
for each set of data that comes in, without having to manually re-type the commands?
Welcome to the beauty and purpose of looping and automating with shell scripts! Read on!

## Loops and bash variables

Right now we've used `grep` and redirecton (`>`) to capture all the bad reads that are present in
one FASTQ file. And you've seen that you can use the grep pattern matching (glob) character
`*` to look at more than one file simultaneously. But using this option works on files
in batch, or all at once, without fine grain, one-by-one control. That's where
looping comes in.

Looping in bash is very similar to other languages. Let's dive right in:

`$ cd ~/unix_oct2015/raw_fastq`

The structure of loops in bash is as follows:

```bash
$ for each in group
> do
>	commands $each
> done
```
where `each` is a variable that takes the value of every member of the specified `group` one at a time and runs through the commands indicated before proceeding to the next member of the group. 

Every loop will have the structure:

```bash
$ for * in *
> do
>   * 
> done
```
You will need to specify the variable name, the group of files, programs, etc. that you would like to perform the command(s) on, and the actual command(s) you would like to perform. For example:

```bash
$ for filename in Mov10_oe_1.subset.fq Mov10_oe_2.subset.fq
> do
>   echo $filename
> done
```

So what does this do? Well, first we specify a group of files: `Mov10_oe_1.subset.fq` and `Mov10_oe_2.subset.fq`. Then, we execute a series of command(s) between the `do` and `done`. But we execute these commands for each item in this list. Moreover, we store in the placeholder, or variable, named `filename` the next FASTQ filename for each time we execute the set of commands. In essence, we loop through the group, executing the commands inside the `do` / `done` block, with the value of `filename` changing each time.

You may have noticed in the `for` statement we reference `filename`. But in the loop, we explicitly use `$filename`. Why? Well, in the former, we're setting the value, while in the latter, we're retrieving the value. This is standard bash notation for setting and getting variables. Forgetting the `$` when you want to retrieve the value of a variable is a common mistake. 

Of course, `filename` is a great variable name. But it doesn't matter what variable name we use:

```bash
$ for x in Mov10_oe_1.subset.fq Mov10_oe_2.subset.fq
> do
>   echo $x
> done
```

The only potential problem is that `x` has little meaning. In the long run, it's best to use a name that will help point out its function, so your future self will understand what you are thinking now.

Looping over two files is great, but its rather inflexible. What notation can we use to grab a whole directory of files? Use the `*` notation!

```bash
$ for filename in *.fq
> do
>   echo $filename
> done
```

Now we've made our looping list more flexible, as it will work on any number of files. So let's put that to work:

```bash
$ for filename in *.fq
> do 
>   echo $filename; 

  # grab all the bad read records
>   grep -B1 -A2 NNNNNNNNNN $filename > ../other/$filename-badreads.fastq
> done
```

In addition to the `echo` statement, we've included a comment -- lines that start with `#` -- and our `grep` statement that finds bad reads and puts them into a new file. And we've specified a new filename with the variable `$filename`. So each iteration of the loop will `grep` a particular file and then output the bad reads to a new file that uses particular filename as part of the new one.

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

