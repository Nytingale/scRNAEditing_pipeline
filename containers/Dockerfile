# scRNAEditing pipeline image: the three conda environments of the benchmark,
# rebuilt from the exact linux-64 lock files in envs/.
#   /opt/conda/envs/scRNAEditing  pipeline, callers and R analyses (on PATH)
#   /opt/conda/envs/scAllele      scAllele 0.0.9.7 (needs numpy 1.26 / scipy 1.15)
#   /opt/conda/envs/cellsnp       Cellsnp-lite 1.2.3 (htslib 1.24)
#   /opt/conda/envs/gt_pileup     samtools 1.20 + bedtools 2.31.0 (ground-truth pileups in analysis/)
#   /opt/conda/envs/refbuild      BEDOPS + bedtools 2.31.0 (extended annotations in references/)
# Build from the repository root:
#   docker build -f containers/Dockerfile -t scrnaediting .
FROM mambaorg/micromamba:2.9.0-debian12-slim@sha256:d8f7ba87cb33f3c35f311748fac5a2c6d28caf807a7507a32dd0cfe1fba13309

USER root
# procps: Nextflow reads task metrics with ps.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*

COPY envs/ /tmp/envs/

RUN micromamba create -y -p /opt/conda/envs/scRNAEditing -f /tmp/envs/scRNAEditing.linux-64.lock.txt \
 && micromamba create -y -p /opt/conda/envs/scAllele -f /tmp/envs/scAllele.linux-64.lock.txt \
 && micromamba create -y -p /opt/conda/envs/cellsnp -f /tmp/envs/cellsnp.linux-64.lock.txt \
 && micromamba create -y -p /opt/conda/envs/gt_pileup -f /tmp/envs/gt_pileup.linux-64.lock.txt \
 && micromamba create -y -p /opt/conda/envs/refbuild -f /tmp/envs/refbuild.linux-64.lock.txt \
 && micromamba clean -a -y

RUN /opt/conda/envs/scAllele/bin/pip install --no-deps --no-cache-dir -r /tmp/envs/scAllele.pip.txt \
 && /opt/conda/envs/scRNAEditing/bin/pip install --no-deps --no-cache-dir setuptools==70.3.0

RUN PATH=/opt/conda/envs/scRNAEditing/bin:$PATH \
    Rscript /tmp/envs/install_r_packages.R /tmp/envs/scRNAEditing.r_packages.tsv \
 && rm -rf /tmp/envs /tmp/Rtmp*

# scAllele and Cellsnp-lite run from their own environments through absolute
# shebangs, so they can sit on PATH next to the main environment.
RUN ln -s /opt/conda/envs/scAllele/bin/scAllele       /usr/local/bin/scAllele \
 && ln -s /opt/conda/envs/scAllele/bin/scAllele_train /usr/local/bin/scAllele_train \
 && ln -s /opt/conda/envs/cellsnp/bin/cellsnp-lite    /usr/local/bin/cellsnp-lite

# What `conda activate scRNAEditing` would set; the JACUSA2 step finds its jar
# through CONDA_PREFIX.
ENV CONDA_PREFIX=/opt/conda/envs/scRNAEditing \
    JAVA_HOME=/opt/conda/envs/scRNAEditing/lib/jvm \
    PATH=/opt/conda/envs/scRNAEditing/bin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

ENTRYPOINT []
CMD ["/bin/bash"]
