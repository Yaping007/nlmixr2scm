# LSF Cheat Sheet

Quick reference for common IBM Spectrum LSF (Load Sharing Facility) commands.

## 1. Checking Host / Cluster Usage

| Command | Description |
|---|---|
| `bhosts` | List all hosts with status and job slot usage (MAX, NJOBS, RUN, etc.) |
| `bhosts -w` | Same as above but with full (wide) status names |
| `bhosts <host>` | Show details for a specific host |
| `lshosts` | Static host info (CPU type, cores, memory, model) |
| `lsload` | Dynamic load info (CPU, memory, paging, I/O) per host |
| `bhosts -l <host>` | Detailed load thresholds and usage for a host |
| `bqueues` | List queues with job slot limits and pending/running counts |
| `bqueues -l <queue>` | Detailed info for a specific queue |

### Reading `bhosts` STATUS
- `ok` — accepting jobs
- `closed` — not accepting new jobs (full, or administratively closed)
- `unavail` / `unreach` — host down or unreachable

### Key columns
- `MAX` — max job slots on the host
- `NJOBS` — total slots used (RUN + SSUSP + USUSP + RSV)
- `RUN` — running jobs
- `SSUSP` / `USUSP` — system- / user-suspended jobs
- `RSV` — reserved slots

## 2. Checking Who Is Using the Cluster

```bash
# Running jobs for all users
bjobs -u all -r

# Per-user count of running jobs (quick summary)
bjobs -u all -r -o "user" | sort | uniq -c | sort -rn

# Custom columns: user, runtime, host, status
bjobs -u all -r -o "user run_time exec_host stat"

# Jobs on a specific host
bjobs -u all -m <host>
```

## 3. Checking Your Job Status

| Command | Description |
|---|---|
| `bjobs` | List your active (pending + running) jobs |
| `bjobs -a` | Include finished jobs |
| `bjobs -r` | Only running jobs |
| `bjobs -p` | Only pending jobs (with reason for pending) |
| `bjobs <jobID>` | Status of a specific job |
| `bjobs -l <jobID>` | Detailed info (resources, host, submit args) |
| `bhist -l <jobID>` | Historical event log for a job |
| `bpeek <jobID>` | View live stdout/stderr of a running job |

## 4. Submitting Jobs

```bash
# Basic submission
bsub -o job.%J.out -e job.%J.err ./myscript.sh

# Request a queue
bsub -q <queue> ./myscript.sh

# Multi-core on a single host
bsub -n 16 -R "span[hosts=1]" -o job.%J.out ./myscript.sh

# Request memory (per slot, MB)
bsub -M 4000 -R "rusage[mem=4000]" ./myscript.sh

# Name a job
bsub -J "myjob" ./myscript.sh

# Interactive session
bsub -Is bash

# Job array (indices 1..100 available via $LSB_JOBINDEX)
bsub -J "myjob[1-100]" -o out.%I.%J.log ./run.sh
```

### Common `bsub` options
- `-n <N>` — number of slots (cores)
- `-R "span[hosts=1]"` — keep all slots on one host
- `-R "span[ptile=<N>]"` — N slots per host
- `-o` / `-e` — stdout / stderr files (`%J` = job ID, `%I` = array index)
- `-q <queue>` — target queue
- `-M <MB>` — memory limit
- `-W <min>` / `-W HH:MM` — wall-clock runtime limit
- `-w "done(<jobID>)"` — dependency (start after another job finishes)

## 5. Managing Jobs

| Command | Description |
|---|---|
| `bkill <jobID>` | Kill a job |
| `bkill 0` | Kill all your jobs |
| `bstop <jobID>` | Suspend a job |
| `bresume <jobID>` | Resume a suspended job |
| `bmod <options> <jobID>` | Modify a pending job's parameters |
| `btop <jobID>` | Move a pending job to the top of the queue |
| `bbot <jobID>` | Move a pending job to the bottom of the queue |

### Changing the concurrency (running-slot limit) of jobs

LSF controls how many jobs/slots run at once via **job limits**. To throttle or
increase how many jobs run concurrently:

```bash
# Job array: limit how many elements run at once with the % modifier.
# Submit an array of 100 that runs at most 10 concurrently:
bsub -J "myjob[1-100]%10" ./run.sh

# Change the concurrency of an ALREADY-SUBMITTED array (e.g. to 20):
bmod -J "myjob[1-100]%20" <jobID>

# Or reference the array by ID only:
bmod "%20" <jobID>        # set running-element limit to 20
bmod "%0"  <jobID>        # 0 = no limit (run as many as slots allow)
```

Notes:
- The `%<n>` after an array name is the **maximum number of running elements** at a time.
- `bmod "%<n>" <jobID>` works only on **job arrays**; it takes effect for
  pending elements immediately (already-running elements keep running).
- Cluster-/queue-/user-wide limits that cap total concurrent jobs are set by
  admins; inspect the ones affecting you with:

```bash
blimits            # show active resource-allocation limits affecting you
blimits -w         # wide output
```

## 6. Cluster & Queue Info

| Command | Description |
|---|---|
| `lsid` | Show cluster name and LSF version |
| `bqueues` | List queues and their limits |
| `busers` | Show per-user job counts and limits |
| `bparams` | Show cluster-wide scheduling parameters |
| `bmgroup` | List host groups |

## 7. Useful One-Liners

```bash
# Count total free slots across ok hosts
bhosts | awk '$2=="ok" {free += $4-$5} END {print free" free slots"}'

# List idle (ok, 0 running) hosts
bhosts | awk '$2=="ok" && $5==0 {print $1}'

# Watch your jobs update every 5s
watch -n 5 bjobs
```
