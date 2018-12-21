
## Run locally
```
export PYSPARK_SUBMIT_ARGS="--driver-memory 8g pyspark-shell"
python 1_run_1000G_test.py
```


## Start dataproc spark server to submit jobs to
```
# Create command using Neale labs cloudtools
cluster start --max-idle 15m --zone europe-west1-b --dry-run em-cluster

# Start server
gcloud beta dataproc clusters create \
    em-cluster \
    --image-version=1.2-deb9 \
    --metadata=MINICONDA_VERSION=4.4.10,JAR=gs://hail-common/builds/0.2/jars/hail-0.2-07b91f4bd378-Spark-2.2.0.jar,ZIP=gs://hail-common/builds/0.2/python/hail-0.2-07b91f4bd378.zip \
    --properties=spark:spark.driver.memory=41g,spark:spark.driver.maxResultSize=0,spark:spark.task.maxFailures=20,spark:spark.kryoserializer.buffer.max=1g,spark:spark.driver.extraJavaOptions=-Xss4M,spark:spark.executor.extraJavaOptions=-Xss4M,hdfs:dfs.replication=1 \
    --initialization-actions=gs://dataproc-initialization-actions/conda/bootstrap-conda.sh,gs://hail-common/cloudtools/init_notebook1.py \
    --master-machine-type=n1-highmem-8 \
    --master-boot-disk-size=100GB \
    --num-master-local-ssds=0 \
    --num-preemptible-workers=0 \
    --num-worker-local-ssds=0 \
    --num-workers=2 \
    --preemptible-worker-boot-disk-size=40GB \
    --worker-boot-disk-size=40 \
    --worker-machine-type=n1-standard-8 \
    --zone=europe-west1-b \
    --initialization-action-timeout=20m \
    --max-idle=15m

# Get lastest hash
HASH=$(gsutil cat gs://hail-common/builds/0.2/latest-hash/cloudtools-3-spark-2.2.0.txt)

gcloud dataproc jobs submit pyspark \
  --cluster=em-cluster \
  --files=gs://hail-common/builds/0.2/jars/hail-0.2-$HASH-Spark-2.2.0.jar \
  --py-files=gs://hail-common/builds/0.2/python/hail-0.2-$HASH.zip \
  --properties="spark.driver.extraClassPath=./hail-0.2-$HASH-Spark-2.2.0.jar,spark.executor.extraClassPath=./hail-0.2-$HASH-Spark-2.2.0.jar" \
  1_run_1000G_dataproc.py

```
