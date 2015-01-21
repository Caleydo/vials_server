caleydo-vis-genome
==================

GenomeViewer (c.web plugin)

## Install

You need to checkout caleydo-web (c.web)

-  `git checkout <url>` into c.web `static/scripts`
-   add to `data/index.json`:

```
    {
        "name": "GenomeBrowser",
        "path": "genomelink.csv",
        "type": "genomeDataLink",
        "idtype":"gdLink"
    },

```

-   add to `data/idtypes.json`:

```
    {
        "id": "gdLink",
        "internal": true,
        "name": "gdlink Columns"
    }
```

- add a file `data/genomelink.csv` with content (2 lines):

```
x,serveradress
0,http://localhost:5000/bam

```

- run `node configgenerator.js` in root of c.web. This will integrate the plugin into
the build file for grunt

- run `grunt build`

- (run the data service provider on port 5000)

- run `grunt serverd` and keep fingers crossed :)

GenomeDataProvider
==================

Reading Genomic Files and provide data as API

run:
python GenomeDataProvider.py


requirements
--------------
manually install matplotlib:
sudo apt-get install python-matplotlib