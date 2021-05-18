# FunctionalMotifs

He creado una carpeta bajo el nombre `Scripts/Scripts_Functional_Groups`. En ella tengo los scripts para calcular los functional groups con el subset de networks muestreadas en un único año (una temporada de floración) y sitio, además son weighted networks con datos de frecuencia de visitas. Con el próposito de los análisis tendría que revisar el resto de redes y añadir las binarias que cumplen esta condición de año y sitio.   

En total son 60 redes de momento y un poco menos de 550 especies de plantas. Los grupos funcionales aparecen de 1 a 5. A ver si añado un gráfico o una tabla que los resuma.  

Todavía estoy trabajando con los scripts así que me queda depurarlos, esto es lo que tengo por ahora:  

*0_Pollinator_Species_Names*: Uso del paquete taxsize para encontrar los nombres de las especies de polinizadores aceptados (los de plantas no están aquí pero seguí el mismo procedimiento).  

*1_Merge_Spp_Names_Long_Format*: En este script remato lo que no encuentra el paquete taxsize de forma manual. En los análisis que estoy haciendo agrego por especie y planta y red los valores de visitas al final del script. _Aquí creo que te interesa los datos brutos sin alterar, así que lo he eliminado y las redes permanecen al completo sin ninguna modificación_.   

*2_Trait_Data_Imputation*: Utilizo dos técnicas de imputación de datos y las dos dan resultados relativamente parecidos (no he visto diferencias en los análisis). Creo que me quedaré con random forest.  

*3_Functional_Groups*: Aquí calculo los grupos funcionales con los datos imputados previamente en el script *2_Trait_Data_Imputation*.  

*4_Merge_Trait_Data_FunctionalG*: Finalmente combino los csv generados en scripts *1_Merge_Spp_Names_Long_Format* y *3_Functional_Groups*. Estos ya serían los datos analizables, quizás habría que quitarle los 0's y singletones. Puedo depurarlos más si quieres y si tienes alguna duda ya me dices! El archivo para analizar está en "Data/Csv/" y se llama "data_for_motifs_analysis.csv".   

Alfonso has extracted motifs `Scripts/Scripts_Functional_Group` and did some analysis `Scripts/Scripts__Functional_Groups_Alfonso`.

- Motifs frequency está en: `Data\Csv\Motifs frequencies and null models\Motifs_frequency_percentile.csv`
- GF_positions está en: `Data\Csv\Motifs positions and null models\GF_positions_frequency_percentile.csv`

Los archivos que he usado para generar los datos de los modelos nulos y procesar las tablas están: `Scripts\Scripts_Alfonso`








