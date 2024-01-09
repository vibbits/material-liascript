<!--
author:   Alexander Botzki
email:    Alexander.Botzki@vib.be
version:  0.2
language: en
narrator: US English Female

comment:  DMP Introduction

icon: https://vib.be/sites/vib.sites.vib.be/files/logo_VIB_noTagline.svg

link:     https://cdnjs.cloudflare.com/ajax/libs/animate.css/3.7.2/animate.min.css
link:     https://raw.githubusercontent.com/vibbits/material-liascript/master/img/org.css
link:     https://elixirtess.github.io/TeSS_widgets/css/tess-widget.css
link:     https://fonts.googleapis.com/css2?family=Saira+Condensed:wght@300&display=swap
link:     https://fonts.googleapis.com/css2?family=Open+Sans&display=swap
link:     vib-styles.css
script:   https://elixirtess.github.io/TeSS_widgets/js/tess-widget-standalone.js" onload="initTeSSWidgets()

debug: true

-->

# scratch

Some text

> ### Questions
>
> Question: A question
>
> > <details markdown="1">
> > <summary>Solution
> > </summary>
> >
> > ### A Header
> > 1. Yes, add explanation here
> >
> > **TODO**: add image
> > </details>

<details markdown="1">
<summary>Solution</summary>
### A Header
1. Yes, add explanation here
**TODO**: add image
</details>

## TeSS Widgets

<div>
<link rel="stylesheet" property="stylesheet" href="https://elixirtess.github.io/TeSS_widgets/css/tess-widget.css"/>
<div id="tess-widget-materials-table" class="tess-widget tess-widget-faceted-table"></div>
<script>
function initTeSSWidgets() {
    TessWidget.Materials(document.getElementById('tess-widget-materials-table'),
        'FacetedTable',
        {
            opts: {
                columns: [{name: 'Name', field: 'title'},
                    {name: 'Description', field: 'description'}],
                allowedFacets: ['scientific-topics', 'target-audience'],
                facetOptionLimit: 5
            },
            params: {
                pageSize: 5,
                q: 'Python'
            }
        });
}
</script>
<script async="" defer="" src="https://elixirtess.github.io/TeSS_widgets/js/tess-widget-standalone.js" onload="initTeSSWidgets()"></script>
</div>
