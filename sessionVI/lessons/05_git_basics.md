---
layout: page
title: Version Control with Git
subtitle: Git basics
minutes: 20
---
## Learning Objectives

*   Understand the benefits of an automated version control system.
*   Explore how to configure Git

## What is version control?

We'll start by exploring how version control can be used
to keep track of what one person did and when.
Even if you aren't collaborating with other people,
automated version control is much better than this situation:

[![Piled Higher and Deeper by Jorge Cham, http://www.phdcomics.com](../img/phd101212s.gif)](http://www.phdcomics.com)

"Piled Higher and Deeper" by Jorge Cham, http://www.phdcomics.com

We've all been in this situation before: it seems ridiculous to have multiple nearly-identical versions of the same document. Some word processors let us deal with this a little better, such as Microsoft Word's "Track Changes" or Google Docs' version history.

Version control systems start with a base version of the document and then save just the changes you made at each step of the way. You can think of it as a tape: if you rewind the tape and start at the base document, then you can play back each change and end up with your latest version.

![Changes are saved sequentially](https://cdn.rawgit.com/hbc/NGS_Data_Analysis_Course/master/sessionVI/img/play-changes.svg)

Once you think of changes as separate from the document itself, you can then think about "playing back" different sets of changes onto the base document and getting different versions of the document. For example, two users can make independent sets of changes based on the same document.

![Different versions can be saved](https://cdn.rawgit.com/hbc/NGS_Data_Analysis_Course/master/sessionVI/img/versions.svg)

If there aren't conflicts, you can even try to play two sets of changes onto the same base document.

![Multiple versions can be merged](https://cdn.rawgit.com/hbc/NGS_Data_Analysis_Course/master/sessionVI/img/merge.svg)

A version control system is a tool that keeps track of these changes for us and
helps us version and merge our files. It allows you to
decide which changes make up the next version, called a
[commit](../reference.html#commit), and keeps useful metadata about them. The
complete history of commits for a particular project and their metadata make up
a [repository](../reference.html#repository). Repositories can be kept in sync
across different computers facilitating collaboration among different people.

> ### The long history of version control systems
>
> Automated version control systems are nothing new. 
> Tools like RCS, CVS, or Subversion have been around since the early 1980s and are used by many large companies.
> However, many of these are now becoming considered as legacy systems due to various limitations in their capabilities.
> In particular, the more modern systems, such as Git and [Mercurial](http://swcarpentry.github.io/hg-novice/) 
> are *distributed*, meaning that they do not need a centralized server to host the repository.
> These modern systems also include powerful merging tools that make it possible for multiple authors to work within 
> the same files concurrently.

## Getting started with Git

When we use Git on a new computer for the first time,
we need to configure a few things. Below are a few examples
of configurations we will set as we get started with Git:

*   our name and email address,
*   to colorize our output,
*   what our preferred text editor is,
*   and that we want to use these settings globally (i.e. for every project)

On a command line, Git commands are written as `git verb`, 
where `verb` is what we actually want to do. So here is how 
Dracula sets up his new laptop:

~~~ {.bash}
$ git config --global user.name "Vlad Dracula"
$ git config --global user.email "vlad@tran.sylvan.ia"
$ git config --global color.ui "auto"
~~~

(Please use your own name and email address instead of Dracula's.)

He also has to set his favorite text editor, following this table:

| Editor             | Configuration command                            |
|:-------------------|:-------------------------------------------------|
| nano               | `$ git config --global core.editor "nano -w"`    |
| Text Wrangler      | `$ git config --global core.editor "edit -w"`    |
| Sublime Text (Mac) | `$ git config --global core.editor "subl -n -w"` |
| Sublime Text (Win, 32-bit install) | `$ git config --global core.editor "'c:/program files (x86)/sublime text 3/sublime_text.exe' -w"` |
| Sublime Text (Win, 64-bit install) | `$ git config --global core.editor "'c:/program files/sublime text 3/sublime_text.exe' -w"` |
| Notepad++ (Win)    | `$ git config --global core.editor "'c:/program files (x86)/Notepad++/notepad++.exe' -multiInst -notabbar -nosession -noPlugin"`|
| Kate (Linux)       | `$ git config --global core.editor "kate"`       |
| Gedit (Linux)      | `$ git config --global core.editor "gedit -s -w"`   |
| emacs              | `$ git config --global core.editor "emacs"`   |
| vim                | `$ git config --global core.editor "vim"`   |

The four commands we just ran above only need to be run once: the flag `--global` tells Git
to use the settings for every project, in your user account, on this computer.

You can check your settings at any time:

~~~ {.bash}
$ git config --list
~~~

You can change your configuration as many times as you want: just use the
same commands to choose another editor or update your email address.

> ### Proxy
>
> In some networks you need to use a
> [proxy](https://en.wikipedia.org/wiki/Proxy_server). If this is the case, you
> may also need to tell Git about the proxy:
>
> ~~~ {.bash}
> $ git config --global http.proxy proxy-url
> $ git config --global https.proxy proxy-url
> ~~~
>
> To disable the proxy, use
>
> ~~~ {.bash}
> $ git config --global --unset http.proxy
> $ git config --global --unset https.proxy
> ~~~
