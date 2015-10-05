---
title: "HPC Introduction"
author: "Bob Freeman, Radhika Khetani"
date: "Monday, October 5, 2015"
---

## Objectives
* Why use a cluster?
* Understand how a cluster is organized
* Understand what a filesystem is and which are appropriate to use
* Know how to interact with the scheduler
* Be comfortable creating a batch script and submitting one
* Know how to get info about jobs and to control them

## TOC
* [Cluster basics](#cluster-basics)
    * [What are some of reasons to access a remote computer system?](#what-are-some-of-reasons-to-access-a-remote-computer-system)
    * [Advantages of using HPC/HTC vs. Cloud systems](#advantages-of-using-hpchtc-vs-cloud-systems)
    * [What does a cluster look like?](#what-does-a-cluster-look-like)
* [Logging in](#logging-in)
* [Filesystems](#filesystems)
* [Using & installing software](#using--installing-software)
* [Working with the scheduler](#working-with-the-scheduler)
    * [Running & submitting jobs](#running--submitting-jobs)
    * [Choosing the proper resources for your job](#choosing-the-proper-resources-for-your-job)
    * [Creating submission scripts](#creating-submission-scripts)
    * [Example batch script (SLURM)](#example-batch-script-SLURM)
    * [Managing jobs and getting job information](#managing-jobs-and-getting-job-information)
* [Best Practices](#best-Practices)
* [Distributed System Definitions and stacks:](#distributed-System-Definitions-and-stacks)
* [HPC vs. Cloud:](#hpc-vs-cloud)
* [Resources:](#resources)

Note: The majority of this document is based on the [Introduction to Odyssey (required) training for HPC users at Harvard University](https://rc.fas.harvard.edu/training/intro-to-odyssey/). Course materials can be found as a PDF file at [https://rc.fas.harvard.edu/wp-content/uploads/2015/05/Intro-to-Odyssey-v2-09.pdf](https://rc.fas.harvard.edu/wp-content/uploads/2015/05/Intro-to-Odyssey-v2-09.pdf)

## Cluster basics

Clusters, otherwise know as high-performance computing (HPC) or high-throughput computing systems; these tools are becoming the <em>de facto</em> standard tools in most research disciplines today.

### What are some of reasons to access a remote computer system?

* Your computer does not have enough resources to run the desired analysis. *E.g.* memory, processors, disk space, or network bandwidth.
* You want to produce results faster than your computer can.
* You cannot install software in your computer. That is, the application does not have support for your operating system, conflicts with other existing applications, or softare licensing does not allow for installation on personal laptops.


### What does a cluster look like?

“High Performance Computing most generally refers to the practice of aggregating computing power in a way that delivers much higher performance than one could get out of a typical desktop computer or workstation in order to solve large problems in science, engineering, or business.” --http://insidehpc.com/hpc-basic-training/what-is-hpc/

Clusters are simply a grouping of computers with the same components (RAM, disk, processors/cores, and networking cards) as those in your desktop or laptop, but with more umph! and are networked with high-speed interconnect that can be accessed (indirectly) through software, the scheduler, that manages simultaneous execution of jobs, or analyses, by multiple persons. 

![Overview of a compute cluster](images/cluster-generic.png)

The user accesses the compute cluster through one or more login nodes, and submits jobs to the scheduler, which will dispatch to and collect the completed work from the compute nodes. Frequently, clusters have shared disks, or filesystems, of various flavors where you can store your data, programs, and use for in-job execution (working or scratch areas)

![Orchestra cluster overview](images/orchestra-outline.png)

## Logging in

ssh
(will ask if you are sure the first time you connect with a remote machine)

ssh -X

scp 


## Filesystems and Storage

What is a filesystem?
Storage on most compute systems is not what and where you think they are! Physical disks are bundled together into a virtual volume; this virtual volume may represent one filesystem, or may be divided up, or partitioned, into multiple filesystems. And your directories then reside within one of these fileystems. Filesystems are accessed over the network through mount points.

![Filesystem definition diagram](images/filesystems-generic.png)
There are often multiple storage/filesystems options available for you to do your work. The most common are:
* home: where you land when you first login
* shared (lab): a directory on the network
* "working" directories, like SCRATCH or TEMP

Home folders are great for keeping information specific to your account, your workflows, and your environment. Typically, these are backed up and are limited in space. Note that these are usually appropriate for small amounts of work -- up to 10 or so jobs -- as these are on filesystems with large #s of other accounts on low-throughput/low bandwidth disk infrastructure.

Shared (lab) systems are typically the same, though may vary from site to site and will vary in size, backup strategy, and usage. These are also usually appropriate for small amounts of work -- up to 10 or so jobs -- as these are on filesystems with large #s of other accounts on low-throughput/low bandwidth disk infrastructure.

"Working directories, often called TEMP, SCRATCH, or WORK, are often specialized, high-availability & high-speed systems designed especially for large volumes of read and write operations found on HPC/HTC systems. Often times this is not backed and files are deleted after an aged period of time -- on Odyssey, this happens for files > 90 days old. Your typical workflow pattern is to stage files in this location prior to or part of your job, do your work, and then copy your final data back to your home or lab share. Ask your local instructor for more guidance on these systems.

Here's a synopsis of filesystems on Orchestra:

!

**Important!! Ensure that you use the proper filesystem for your work, as improper usage can negatively affect other people’s jobs on those same physical disks.**

**Exercises**
* 

## Using & installing software

On most cluster systems, all the software installed is not immediately available for use;
instead, it is loaded into your environment incrementally using a module system. And our module system is pretty easy to use:

<pre>
$ module load seq/fastqc				# get the lastest version of a program
$ module list							# list the modules that have been loaded into your environment
$ module load seq/fastqc/0.10.1		# load a specific version
</pre>

To find software installed and available on the cluster, two commands can help you:

<pre>
$ module avail						# list all available modules
$ module avail seq/					# list all available modules related to sequence analysis
$ module avail seq/fastq			# list all available modules for tools starting with the word "fastq"
</pre>

|| *For Perl & Python modules or R packages*, we encourage you to set up directories in your
home and/or lab folders for installing your own copies locally. Please see [[our handy
instructions]](instructions available?) for more information.

|| *If software you need is not installed*, we encourage you to do local installs in your home
or lab folder for bleeding-edge releases, software you are testing, or software used
only by your lab. See [[our helpful guide]](link to guide) for more information. For programs that are commonly used by your domain, field, or department, please submit a 
[[software install request]](request submission link).

|| Note that due to demand and the complex nature of software installs, it may take one to two 
weeks for us to complete these requests.


## Working with the scheduler

As mentioned before, the scheduler is responsible for listening to your job requests, then finding the proper compute node that meets your job's resource requirements -- RAM, # cores, time, etc -- dispatches the job to that compute node, collects info about the completed work, and stores information about your job. If you've asked it to do so, it will even notify you about the status of your job (e.g. begin, end, fail, etc).


### Running & submitting jobs

There are two ways to run jobs on a cluster. One, usually done at the start, is to get an
interactive session on a compute node. This looks and behaves exactly as when you first log into a compute
cluster, on a login node, but the work is being done a compute node, a worker node on the cluster. This is 
considered a best practice technique, and should be done for all work that will tie up resources (e.g. CPU-
or memory-intensive tasks).

To get an interactive session, you issue the `bsub` command with the appropriate parameters for requesting
the resources you require. For example:

```bash
$ bsub -Is -q interactive bash
```

This command requests from the scheduler a foreground/interactive job with the following resources:
```
-Is				# submits a batch interactive job and creates a pseudo-terminal with 
					#shell mode support when the job starts
-q interactive	# the queue to run on
bash			# the program we want to run, which is the bash shell
```

A few additional, optional, parameters were left out; as such, LSF will give us the defaults:
```bash
-n 1            		# number of cores (CPUs), default = 1
-R "rusage[mem=2000]"	# amount of memory, default = 2GB
[[]]-W 						# amount of "wall clock time", default = ??
```

The other way is to create a batch submission script file, which has these parameters embedded inside, and submit your script to the scheduler:

```bash
bsub < my_batch_script.sh
```

In all cases, the scheduler will return to you a jobID, a unique ID for your job that you can use to get info or control at that time, or refer to it historically.

### Choosing the proper resources for your job
For both types of submissions, you are requesting resources from the scheduler to run your job. These are:
* total time (wall clock time)
* memory
* # of cores (CPUs)
* # of nodes
* a queue 

Choosing resources is like playing a game with the scheduler: You want to request enough to get your job completed without failure, But request too much: your job is ‘bigger’ and thus harder to schedule. Request too little: if your job goes over that requested, it is killed. So you want to get it just right, and pad a little for wiggle room.

Another way to think of 'reserving' a compute node for you job is like making a reservation at a restaurant:
* if you bring more guests to your dinner, there won't be room at the restaurant, but the wait staff may try to fit them in. If so, things will be more crowded and the kitchen (and the whole restaurant) may slow down dramatically.
* if you bring fewer guests and don't notify the staff in advance, the extra seats are wasted; no one else can take the empty places, and the restaurant may lose money.

“Never use a piece of bioinformatics software for the first time without looking to see what command-line options are available and what default parameters are being used”
	-- acgt.me · by Keith Bradnam
	
#### Time
This is determined by test runs that you do on your code during an interactive session.
Or, if you submit a batch job, over-ask first, check the amount of time actually needed,
then reduce time on later runs.


#### Memory:
|| We recommend that you check the software docs for memory requirements. But often times these are not stated, so we can take another approach. On Orchestra, each job is allowed, on average, 4 GB RAM/core allocated. So, try 4 GB and do a trial run. If your job was killed, look at your log files or immediately with squeue (see later). If it shows a memory error, you went over. Ask for more and try again.

|| Once the job has finished, ask the scheduler how much RAM was used by using the `sacct` command to get post-run job info:
```bash
sacct -j JOBID --format=JobID,JobName,ReqMem,MaxRSS,Elapsed	# RAM requested/used!!
```
|| The `ReqMem` field is how much you asked for and `MaxRSS` is how much was actually used. Now go back and adjust your RAM request in your sbatch command or submission script.


#### # of Cores
This is determined by your software, how anxious you are to get the work done, and how well your code scales. **NOTE! Throwing more cores at a job does not always make it run faster!** This is often a newbie mistake and will waste compute, making your admin grumpy. Ensure your software can use multiple cores: Inspect the parameters for your software and look for options such as 'threads', 'processes', 'cpus', 'cores'; this will often indicate that it has been parallelized. Then run test jobs to see how well it performs with multiple cores, inching slowing from 1 to 2, 4, 8, etc, assessing the decrease in time for the job run as you increase cores. Parallelizable programs often do not scale well -- it's important to understand this so you can choose the appropriate number.

#### # of Nodes
For most software in biology, this choice is simple: 1. There are *very* few biology softare packages capable of running across multiple nodes. If they are capable, they will mention the use of technology called 'MPI' or 'openMPI'. Please talk to your local, friendly Research Computing Facilitator (e.g. XSEDE Campus Champion) on how to make use of this feature.

#### Queues

Queues, are a grouping of computers to run a certain profile of jobs. This could be maximum run time, # of cores used, amount of max RAM, etc. For Orchestra there is a [handy guide](https://wiki.med.harvard.edu/Orchestra/ChoosingAQueue) to choosing the correct queue. Queues will vary significantly from site to site, so please check with your instructor for the appropriate ones to use, and why.

Here's a synopsis of partitions unique to Orchestra:

### Creating submission scripts

When creating submission scripts, use your favorite text editor and include the following sections:

**Required**
* shebang, or shell invocation
* Partition to submit to (comma separated)
* Any sofware or module loads
* actual commands to do work

**Recommended**
* job name
* # of cores
* # of nodes
* amount of time
* amount of memory
* STDOUT file
* STDERR file

**Optional**
* what mail notifications you want, if any
* email notification address (or text msg email address if you want an SMS alert!)

All the scheduler directives need to be at the start of the file, then your module/software loads, and then your actual job commands. 

### Example batch script (SLURM)

The following is an example submission script for an example `fastqc` run, which is a single-core program:

```bash
#!/bin/bash
#
#SBATCH -p serial_requeue        # Partition to submit to (comma separated)
#SBATCH -J frog_fastqc           # Job name
#SBATCH -n 1                     # Number of cores
#SBATCH -N 1                     # Ensure that all cores are on one machine
#SBATCH -t 0-0:10                # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 1000               # Memory in MB
#SBATCH -o fastqc_%j.out         # File for STDOUT (with jobid = %j)
#SBATCH -e fastqc_%j.err         # File for STDERR (with jobid = %j)
#SBATCH --mail-type=ALL          # Type of email notification: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=EMAIL@U.COM  # Email where notifications will be sent

source new-modules.sh; module load fastqc

cd my_output_directory

fastqc --casava -o fastqc_reports A1_R1.pair.fastq.gz

#... more processing goes here ...

```
Please refer to the following for other examples of SLURM scripts:
* [multicore job script](https://rc.fas.harvard.edu/wp-content/uploads/2015/05/Intro-to-Odyssey-v2-09.pdf#page=34)
* [multicore openMP job script](https://rc.fas.harvard.edu/wp-content/uploads/2015/05/Intro-to-Odyssey-v2-09.pdf#page=35)
* [multinode/MPI job script](https://rc.fas.harvard.edu/wp-content/uploads/2015/05/Intro-to-Odyssey-v2-09.pdf#page=36)

**Exercises**
* What would you need to change to make this a multicore (e.g. 8) job? Hint: more than one change is needed!
* Create a SLURM script from scratch, asking to run `bwa mem` on all your trimmed FASTQs on the general partition, with 2 GB of RAM, 1 core, for 1 hr.
* Create another SLURM script, now asking to run this same script using 4 cores and the parallel option for `bwa mem`.
* (Bonus!) If you had completed the script in the [Shell: Loops & Scripts](https://github.com/fasrc/DataC-shell-genomics/blob/gh-pages/lessons/03_loops_and_scripts.md) lesson, create a new script based on  your `generate_bad_reads_summary.sh` that contains the appropriate SLURM directives. When finished, submit this job to the cluster.

### Managing jobs and getting job information

There are several commands that you can use to control and get info about your jobs:

`scancel` will become your friend! At some point, you'll fire off one or more jobs, and realize you've
made a mistake. (What? You don't make them? Then you can forget about this command) Here are a 
few examples of `scancel` in action:

```bash
scancel JOBID										# specific job
scancel -u bfreeman								    # ALL my jobs
scancel -u bfreeman -J many_blast_jobs				# named jobs
scancel -u bfreeman -p bigmem						# ALL in partition
```
`squeue` will give you pending (to be done), running, and recently completed job info. Some examples:

```bash
squeue -u bfreeman									# jobs for bfreeman
squeue -u bfreeman --states=R | wc –l				# # of Running jobs
```

`sacct` will give you current and historical information, since time began or you were an HPC-infant,
whichever came first. More examples:

```bash
sacct -u bfreeman									# jobs for bfreeman
sacct -u bfreeman -p bigmem --starttime=9/1/14		# same+bigmem partition
sacct -j JOBID --format=JobID,JobName,ReqMem,MaxRSS,Elapsed # RAM requested & used!!
```

**Other info**
See our list of common SLURM commands at https://rc.fas.harvard.edu/resources/documentation/convenient-slurm-commands/

For those coming from another cluster/scheduler, check our our scheduler Rosetta stone: http://slurm.schedmd.com/rosetta.pdf

**Exercises**
At the start of these exercises, your instructor will start a number of jobs. Once started, please
complete the following questions:
* What `scancel` command would you issue to cancel pending jobs in the serial_requeue partition?
* What `squeue` command would you use to count how many jobs are in each state?
* What `sacct` command would you use to see all your jobs named BLAST run between Aug 1st and Aug 15th?
* (Bonus!) Give the command to ....


## Best Practices

Again, working on a cluster is working in a big sandbox, with people of all ages and skills. So it is
important to work carefully and be considerate. Please visit our list of Common Pitfalls and 
Fair Use/Responsibilities pages so that you'll be a good member of the community...

Common Pitfalls: https://rc.fas.harvard.edu/resources/documentation/common-odyssey-pitfalls/<br>
Fair Use/Responsibilities: https://rc.fas.harvard.edu/resources/responsibilities/


## Distributed System Definitions and stacks:
  (Note that many definitions exist for these terms)

 * Distributed application: an application that can be executed on a distributed system platform (e.g., mpiBLAST)
 * Distributed system platform: software layers that facilitates coordination and management of a distributed system (e.g., queue-based system, and MapReduce)
 * Distributed system:
   * High Performance Computing (HPC): large assemble of physical machines and a homogeneous operating system (e.g., your institutions' HPC, XSEDE's HPC)
   * Cloud Computing: virtual machines, distributed platforms and/or applications offered as a service (e.g., Amazon Web Services, Microsoft Azure, Google Cloud Computing)

 * Virtual machine (VM): software computer that like a physical computer can run an operating system and applications
 * Operating system (OS): the basic software layer that allows execution and management of applications
 * Physical machine: the hardware (processors, memory, disk and network)

## HPC vs. Cloud:

| HPC | Cloud |
|:----|:------|
| User account on the system | root account on the system |
| Limited control of the system | Full control of the system |
| Central shared file system | Local file system |
| Jobs submitted into a queue | Jobs executed on each resource |
| Account-based isolation | OS-based isolation |
| Batch-oriented execution of applications | support for batch or interactive applications |
| Request for resource and time allocation | Pay-per-use |
| etc. | etc.|

![HPC vs. Cloud](https://raw.githubusercontent.com/datacarpentry/cloud-genomics/master/lessons/images/HpcVsCloud.png)

### Advantages of using HPC/HTC vs. Cloud systems

* no need to transfer files before and after a "remote" session
* security in knowing that you control how your data is managed
* access to local expertise
* economies of scale
* cost: often local resources may be free or you'll know the cost up-front
* cloud work is well-suited for highly-parallelized, independent workloads (non-MPI)

## Resources:
 * HPC offerings:
  * Odyssey: https://rc.fas.harvard.edu/odyssey/
  * XSEDE: https://www.xsede.org/high-performance-computing
 * Cloud computing offerings:
  * Amazon EC2: http://aws.amazon.com/ec2/
  * Microsoft Azure: https://azure.microsoft.com/en-us/
  * Google Cloud Platform: https://cloud.google.com/
  * iPlant's Atmosphere: http://www.iplantcollaborative.org/ci/atmosphere
