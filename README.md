<img src="./docs/pics/DEUS_logo_v2.png" alt="DEUS logo" width="300px" /> 

# An R package for accurate small RNA profiling based on differential expression of unique sequences 

Tim Jeske, Peter Huypens, Laura Stirm, Selina Höckele, Christine Wurmser, Anja Böhm, Cora Weigert, Harald Staiger, Christoph Klein, Johannes Beckers, and Maximilian Hastreiter

The DEUS R package is a framework for accurate small non-coding RNA (sncRNA) profiling based on differential expression of unique sequences. In comparison to commonly-used 
mapping-based approaches, DEUS generates read counts on the actual sequence reads. Several advantages of this approach include:

- Reads that do not map or map to multiple features are retained in the analysis and represented in an unprecedented and more concise manner
- Sequence clustering provides unique insight into potential processing and editing steps among members of the same cluster
- No need for a reference genome and, therefore, facilitates sncRNA profiling in virtually any organism
- Easy to use and flexible due to high modularity
- Returns intuitively interpretable results

Full documentation is available at <https://timjeske.github.io/DEUS/>

To directly use DEUS you can download and run our <a href="https://github.com/timjeske/DEUS/raw/master/sample_usage.R"> test script</a> which comes along with its own test data.
