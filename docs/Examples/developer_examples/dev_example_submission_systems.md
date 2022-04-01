---
date: '2022-04-01T13:55:18.499Z'
docname: Examples/developer_examples/dev_example_submission_systems
images: {}
path: /examples-developer-examples-dev-example-submission-systems
title: Submission Systems
---

# Submission Systems

Submission system play an important role, if you want to develop your pygromos code. Many times, they are hidden in the Simulation_runner blocks. But maybe you want to develop something, where you need direct access on the submission system?

This notebook will give you some examples, how you can use the submission systems. Note that all submission systems are write in the same ways, such you can exchange them quickly.

## Local Submission

This system executes the commands directly in your current session. This allows you to locally test or execute your code. Maybe if your process needs much more time, you want later to switch to a submission system for job-queueing.

## LSF Submission

The Lsf submission system allows to submit jobs to the IBM LSF-Queueing system.

**Careful! This part requires a running LSF-Queueing System on your System**

You can submit and kill jobs and arrays to the queue, as well as getting information from the queuing list.

### Queue Checking:

### Submission:

here you can submit jobs to the queue as bash commands

### Submitting multiple jobs

### Killing a jobs

Remove a job the job queue
