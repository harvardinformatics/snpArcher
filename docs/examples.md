# Examples
On this page you will find example project scenarios and how to setup and run them using snpArcher.

In this example, we have 10 resequenced individuals we would like to generate variant calls for. We will cover creating the sample sheet, selecting config options, and running the workflow.
## Directory structure
First, let's setup our directories as suggested in our [executing](./executing.md#optional-directory-setup) instructions. Let's assume we are working in a directory called `workdir/`, and the snpArcher repository has already been cloned there. We have also already created the `snparcher` conda env as instructed in the [setup docs](./setup.md#environment-setup).

1. Let's create a directory to organize this project and future ones, call it `projects`. Then, create a new directory for this project, we'll call it `secretarybird_reseq`. 
```
.
├── projects
│   └── secretarybird_reseq
└── snpArcher
```
```{note}
Not all files and directories are shown, only relevant ones. 
```
2. Copy the snpArcher config directory `snpArcher/config` to `projects/secretarybird_reseq`:
```
.
├── projects
│   └── secretarybird_reseq
│       └── config
│           └── config.yaml
└── snpArcher
```

3. Assume we already have all our sequence data and reference genome on our system, stored in a different location `/storage/data`. We do not need to move the raw data to our project directory. 
```{note}
We'll cover the cases using SRA data and refSeq genomes later on in this example.
```