---
Category: literaturenote
Year: {{date | format("YYYY")}}
Tags: literature {% if allTags %}{{allTags}}{% endif %} 
Citekey: {{citekey}}
Authors: {{authors}}{{directors}}
---
[[References]] 

>[!Cite] 
> {{bibliography}} 
  
>[!synth] 
>**Contribution**::  
>  
>**Related**:: {% for relation in relations | selectattr("citekey") %} [[@{{relation.citekey}}]]{% if not loop.last %}, {% endif%} {% endfor %} 
>

  
>[!md]  
{% for type, creators in creators | groupby("creatorType") -%}  
{%- for creator in creators -%}  
> **{{"First" if loop.first}}{{type | capitalize}}**:: 
{%- if creator.name %} {{creator.name}}  
{%- else %} {{creator.lastName}}, {{creator.firstName}}  
{%- endif %}  
{% endfor %}~  
{%- endfor %}  
> **Title**:: {{title}} 
> **Year**:: {{date | format("YYYY")}} 
> **Citekey**:: {{citekey}} {%- if itemType %} 
>**itemType**:: {{itemType}}{%- endif %}{%- if itemType == "journalArticle" %}  
> **Journal**::*{{publicationTitle}}*== {%- endif %}{%- if volume %}  
>**Volume**:: {{volume}} {%- endif %}{%- if issue %} 
>**Issue**:: {{issue}} {%- endif %}{%- if itemType == "bookSection" %}  
>**Book**:: {{publicationTitle}} {%- endif %}{%- if publisher %}  
> **Publisher**:: {{publisher}} {%- endif %}{%- if place %} 
>**Location**:: {{place}} {%- endif %}{%- if pages %}  
>**Pages**:: {{pages}} {%- endif %}{%- if DOI %}  
>**DOI**:: {{DOI}} {%- endif %}{%- if ISBN %}  
> **ISBN**:: {{ISBN}} {%- endif %} 
  

  
> [!Abstract]  
> {%- if abstractNote %} 
> {{abstractNote}}  
> {%- endif -%}.
>


## Annotations
{% for annotation in annotations -%} 
    {%- if annotation.annotatedText -%} 
    {{annotation.annotatedText}}”{% if annotation.color %} {{annotation.colorCategory}} {{annotation.type | capitalize}} {% else %} {{annotation.type | capitalize}} {% endif %}[Page {{annotation.page}}](zotero://open-pdf/library/items/{{annotation.attachment.itemKey}}?page={{annotation.page}}&annotation={{annotation.id}}) 
    {%- endif %} 
    {%- if annotation.imageRelativePath -%}
    ![[{{annotation.imageRelativePath}}]] {%- endif %} 
{% if annotation.comment %} 
{{annotation.comment}} 
{% endif %} 
{% endfor -%}

--- 

## Notes
See [[{{citekey}}-notes]] for manual notes.

















