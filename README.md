# GenonaBuilding
Bio-informatic assignment, we done a trimming operation on FASTQ file, built the overlap graph and de Bruijn graph for genona rebuilding. We can manage bubbles and tips error on de Bruijn graph, we assume that we see what is correct several times, compared to what is not; we also compared the two graph on memory usage and time request.

The documentation and the rest of README is in italian for our convenience; the code is english commented in its entirety. 

<p align="center">
  <img src="https://user-images.githubusercontent.com/24355671/93601343-59463880-f9c1-11ea-81df-bb05bb9d79ab.gif">
</p>

# Assemblaggio di reads

Il progetto consiste in una pipeline di assemblaggio di read provenienti da un genoma, che dato in input un file in formato FASTQ,
esegue tre step in modo sequenziale. Il primo step è il trimmming dell'insieme di reads, secondo step consiste nella costruzione di un grafo di Overlap e un grafo di De Bruijn e come terzo step si effettua la visita dei grafi per produrre in output il genoma finale.
### Prerequisites

Strumenti necessari per il corretto funzionamento.

```
Python 3.7.3
Python3-pip
datetime - python library
difflib  - python library
memory_profiler - python library
numpy - python library 
igraph - python library
operator - python library
```
### Installing

Steps per l'installazione degli strumenti.

```
pip install DateTime
pip install memory-profiler
pip install numpy
pip install igraph
pip install operator-courier
```
## Running the code
All'interno della cartella Assignment eseguire :
```
python ./Code/Main.py 
```

### Break down into end to end tests

```
I risultati si possono verificare in ProgettoBioinformtica/Assignment/File:
CleanReadsNGS.txt
FinalGenomaByDeBrujin.txt
FinalGenomaByOverlap.txt
```

## Deployment

```
git clone https://github.com/BrunoPalazzi42/ProgettoBioinformtica.git
```

## Built With

* [Python 3.7.3](https://www.python.org/downloads/) - linguaggio utilizzato.

## Authors

* **Villa Giacomo** - *g.villa48@campus.unimib.it* - [Villone96](https://github.com/Villone96)
* **Palazzi Bruno** - *b.palazzi@campus.unimib.it* - [BrunoPalazzi42](https://github.com/BrunoPalazzi42)
* **Matamoros Ricardo** - *r.matamorosaragon@campus.unimib.it* - [ricardo_aragon](https://github.com/ricardoanibalmatamorosaragon)
