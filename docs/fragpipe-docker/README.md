# FragPipe Docker Containers for nf-core

## ⚠️ NO REDISTRIBUTION — READ THIS FIRST ⚠️

# **You MUST NOT redistribute FragPipe Docker images.**

# **FragPipe tools are licensed software. Uploading built images to Docker Hub, GitHub Container Registry, Quay.io, or any other public registry is a violation of the FragPipe license agreement.**

**After building, store your image in:**
- A **local Docker daemon** (default after `docker build`)
- A **private container registry** (e.g., AWS ECR, private Docker Hub repo, Artifactory)

---

## Overview

nf-fragpipe supports two Docker container variants:

- **Academic container** — For academic and non-profit use (free license from NESViLab)
- **Commercial container** — For industry use (paid license from Fragmatics)

Both containers use `fcyucn/fragpipe:24.0` as their base image, which provides all open-source FragPipe tools pre-installed (Philosopher, Percolator, MSBooster, PTMShepherd, etc.).
The academic container layers academic-licensed JARs (MSFragger, IonQuant, DiaTracer) on top.
The commercial container overlays the FragPipePlus distribution (commercial JARs + DIA-NN) on the same base.

Both containers expose the same tool aliases in PATH, so pipeline modules work identically with either variant.

## Academic Container (Step-by-Step)

### Prerequisites

- Docker installed and running
- Internet access to download tools

### Step 1: Download Licensed Tools

You need three tools from NESViLab (free for academic/non-profit use):

1. **MSFragger**: <https://msfragger-upgrader.nesvilab.org/upgrader/>
   - Fill in your academic email and agree to the license
   - Download the zip (e.g., `MSFragger-4.4.1.zip`)

2. **IonQuant**: <https://msfragger-upgrader.nesvilab.org/ionquant/>
   - Fill in your academic email and agree to the license
   - Download the zip (e.g., `IonQuant-1.11.20.zip`)

3. **DiaTracer**: <https://msfragger-upgrader.nesvilab.org/diatracer/>
   - Fill in your academic email and agree to the license
   - Download the zip (e.g., `diatracer-2.2.1.zip`)

### Step 2: Build the Container

```bash
cd nf-fragpipe/docker

export MSFRAGGER_ZIP=~/Downloads/MSFragger-4.4.1.zip
export IONQUANT_ZIP=~/Downloads/IonQuant-1.11.20.zip
export DIATRACER_ZIP=~/Downloads/diatracer-2.2.1.zip

./build_nfcore.sh academic
```

### Step 3: Verify

```bash
docker run --rm fragpipe-nfcore:24.0 msfragger --version
docker run --rm fragpipe-nfcore:24.0 ionquant --version
docker run --rm fragpipe-nfcore:24.0 diatracer --version
docker run --rm fragpipe-nfcore:24.0 philosopher version
docker run --rm fragpipe-nfcore:24.0 percolator --help 2>&1 | head -1
```

### Step 4: Save to Private Registry

```bash
# Tag for your private registry
docker tag fragpipe-nfcore:24.0 your-registry.example.com/fragpipe:24.0-academic

# Push to your private registry
docker push your-registry.example.com/fragpipe:24.0-academic
```

### Step 5: Use with Pipeline

```bash
nextflow run main.nf \
  --mode fragpipe \
  --fragpipe_container your-registry.example.com/fragpipe:24.0-academic \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -profile docker
```

## Commercial Container (Step-by-Step)

### Prerequisites

- Docker installed and running
- Commercial license from Fragmatics

### Step 1: Obtain Commercial License

1. Contact **Fragmatics** at info@fragmatics.com or visit <https://www.fragmatics.com/>
2. Purchase a commercial license for FragPipePlus
3. You will receive:
   - `FragPipePlus-24.0-linux.zip` — the complete distribution
   - `license.dat` — the license file (required at runtime by MSFragger, IonQuant, DiaTracer)

### Step 2: Build the Container

```bash
cd nf-fragpipe/docker

export FRAGPIPEPLUS_ZIP=~/Downloads/FragPipePlus-24.0-linux.zip
export LICENSE_FILE=~/Downloads/license.dat

./build_nfcore.sh commercial
```

### Step 3: Verify, Save & Use

Same as academic container (Steps 3-5 above). Replace `academic` with `commercial` in the tag.

## Custom Image Tag

Override the default tag with the `TAG` environment variable:

```bash
TAG=my-fragpipe:latest ./build_nfcore.sh academic
```

## Container Override

The root pipeline defaults to `fcyucn/fragpipe:24.0` (open-source base, no licensed tools).
Use `--fragpipe_container` to switch to your locally built image:

```bash
nextflow run main.nf \
  --mode fragpipe \
  --fragpipe_container your-registry.example.com/fragpipe:24.0-academic \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -profile docker
```

All FragPipe-based workflows (DDA LFQ, TMT Label Check, generic FragPipe) respect this parameter.

## Tool Versions

| Tool | Version | License |
|------|---------|---------|
| Java | 17 (openjdk-17, from base image) | GPL-2.0 |
| MSFragger | 4.4.1 | Academic/Commercial |
| IonQuant | 1.11.20 | Academic/Commercial |
| DiaTracer | 2.2.1 | Academic/Commercial |
| Philosopher | 5.1.3-RC9 | GPL-3.0 |
| Percolator | 3.7.1 | Apache-2.0 |
| MSBooster | 1.4.14 (academic) / 1.4.17 (commercial) | LGPL-3.0 |
| PTMShepherd | 3.0.11 | Open Source |
| TMT-Integrator | 6.1.3 | Open Source |
| Crystal-C | 1.5.10 | Open Source |
| DIA-Umpire | 2.3.3 | Open Source |
| MBG | 0.3.6 | Open Source |
| batmass-io | 1.36.5 (academic) / 1.36.6 (commercial) | Open Source |
| DIA-NN | 1.8.1 | Free (bundled in FragPipe container) |

## Architecture

Both containers are based on `fcyucn/fragpipe:24.0`, which installs the full FragPipe distribution at `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/`.
Pipeline modules discover tools at runtime via the `tools_dir` path (default: `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools`).

```
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/
├── bin/                  (FragPipe CLI)
├── lib/                  (FragPipe Java libraries)
├── tools/
│   ├── MSFragger*.jar    (licensed, added by academic/commercial build)
│   ├── IonQuant*.jar     (licensed, added by academic/commercial build)
│   ├── DiaTracer*.jar    (licensed, added by academic/commercial build)
│   ├── MSBooster*.jar
│   ├── Philosopher/philosopher  (binary)
│   ├── percolator        (binary)
│   ├── PTMShepherd*.jar
│   ├── TMTIntegrator*.jar
│   ├── CrystalC*.jar
│   ├── DIAUmpire*.jar
│   ├── MBG*.jar
│   └── ...
└── ext/
    ├── bruker/           (Bruker native libs)
    └── thermo/           (Thermo native libs)
```

## Troubleshooting

### `msfragger: command not found`

Ensure the container was built correctly and the licensed JAR is in the tools directory.
Check that the JAR exists:

```bash
docker run --rm <image> find /fragpipe_bin -name "MSFragger*.jar"
```

### `No MSFragger*.jar found`

The licensed JAR was not properly extracted or copied into the tools directory.
For academic builds, verify that your `MSFRAGGER_ZIP` points to a valid zip containing the JAR.
Check that the zip structure matches (e.g., `MSFragger-4.4.1/MSFragger-4.4.1.jar`):

```bash
unzip -l $MSFRAGGER_ZIP | head -20
```

### JAVA_HOME mismatch

Both containers inherit Java 17 (`openjdk-17`) from the `fcyucn/fragpipe:24.0` base image.

If you see errors like `UnsupportedClassVersionError` or `JAVA_HOME is not set`, check which Java is active:

```bash
docker run --rm <image> java -version
docker run --rm <image> bash -c 'echo $JAVA_HOME'
```

FragPipe 24.0 tools are compatible with Java 17+.
Do not downgrade to Java 11 or earlier, as MSFragger 4.x requires Java 17 at minimum.
