# MetaboKit pipeline
## DDA quantitation and library generation
```
cd {path_to_folder_with_mzML_files_and_param.txt}
python3 {path_to_metabokit}/linux/DDAfeature.py
python3 {path_to_metabokit}/linux/DDAscore.py
python3 {path_to_metabokit}/linux/DDAalign.py
```

## DIA quantitation
```
python3 {path_to_metabokit}/linux/DIAfeature.py
python3 {path_to_metabokit}/linux/DIAscore.py
python3 {path_to_metabokit}/linux/DIAalign.py
```


