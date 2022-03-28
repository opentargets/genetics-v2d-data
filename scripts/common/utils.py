from pathlib import Path
from psutil import virtual_memory

from pyspark.conf import SparkConf
from pyspark.sql import SparkSession


def get_project_root() -> Path:
    return Path.cwd()


def get_config_path() -> Path:
    return 'configs' / 'config.yaml'


def detect_spark_memory_limit():
    """Spark does not automatically use all available memory on a machine. When working on large datasets, this may
    cause Java heap space errors, even though there is plenty of RAM available. To fix this, we detect the total amount
    of physical memory and allow Spark to use (almost) all of it."""
    mem_gib = virtual_memory().total >> 30
    return int(mem_gib * 0.9)


def initialize_sparksession() -> SparkSession:
    """Initialize spark session."""

    spark_mem_limit = detect_spark_memory_limit()
    spark_conf = (
        SparkConf()
        .set('spark.driver.memory', f'{spark_mem_limit}g')
        .set('spark.executor.memory', f'{spark_mem_limit}g')
        .set('spark.driver.maxResultSize', '0')
        .set('spark.debug.maxToStringFields', '2000')
        .set('spark.sql.execution.arrow.maxRecordsPerBatch', '500000')
    )
    spark = (
        SparkSession.builder.config(conf=spark_conf)
        .master('local[*]')
        .config("spark.driver.bindAddress", "127.0.0.1")
        .getOrCreate()
    )

    return spark
