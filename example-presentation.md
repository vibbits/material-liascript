<!--
title: "Liascript Presentations"

import: https://raw.githubusercontent.com/LiaScript/CodeRunner/master/README.md

icon: https://tess.elixir-europe.org/assets/elixir/elixir-tess-219b707c4912e9c46c917a24ce72b464ec9f2fd56ce03dbcee8b2f6b9ac98a44.svg

link:   https://cdnjs.cloudflare.com/ajax/libs/animate.css/4.1.1/animate.min.css
link: https://fonts.googleapis.com/css?family=Lato:400,400italic,700
link: style.css

@runR: @LIA.eval(`["main.R"]`, `none`, `Rscript main.R`)
-->

# Hello, There

This presentation will show you examples of what you can do with [Liascript](https://liascript.github.io), including:

-   Presenting code and LaTeX equations
-   Including computations in slide output
-   Image, video, and iframe embedding
-   Elegant transitions and animations
-   Printing to PDF via Liascript Exporter

...and much more

## Pretty Code

-   Many syntax highlighting themes available[^1]
-   Three subsequent backticks indicate a stack code block

Example for an R script

``` R
# Define a server for the Shiny app
function(input, output) {
  
  # Fill in the spot we created for a plot
  output$phonePlot <- renderPlot({
    # Render a barplot
  })
}
```

[^1]: Learn more about [Syntax Highlighting](https://github.com/LiaScript/docs/blob/master/Code.md)

## Code Blocks

- grouped code-blocks are simply attached to each other[^1]
- code-blocks can be visible or hidden

``` js     -EvalScript.js
let who = data.first_name + " " + data.last_name;

if(data.online) {
  who + " is online"; }
else {
  who + " is NOT online"; }
```
``` json    +Data.json
{
  "first_name" :  "Sammy",
  "last_name"  :  "Shark",
  "online"     :  true
}
```

[^1]: Learn more about [Projects](https://liascript.github.io/course/?https://raw.githubusercontent.com/liaScript/docs/master/README.md#41)

## Executable Code

``` R
#| echo: true
#| fig-width: 10
#| fig-height: 4.5
library(ggplot2)
png(file="out2.png")

ggplot(mtcars, aes(hp, mpg, color = am)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "loess")

```
@runR

[^1]: Learn more about [Interactive code](https://liascript.github.io/course/?https://raw.githubusercontent.com/liaScript/docs/master/README.md#41)

## LaTeX Equations

[MathJax](https://www.mathjax.org/) rendering of equations to HTML[^1]


$$
\begin{gather*}
a_1=b_1+c_1\\
a_2=b_2+c_2-d_2+e_2
\end{gather*}

\begin{align}
a_{11}& =b_{11}&
  a_{12}& =b_{12}\\
a_{21}& =b_{21}&
  a_{22}& =b_{22}+c_{22}
\end{align}
$$

$$
\begin{gather*}
a_1=b_1+c_1\\
a_2=b_2+c_2-d_2+e_2
\end{gather*}
$$
$$
\begin{align}
a_{11}& =b_{11}&
  a_{12}& =b_{12}\\
a_{21}& =b_{21}&
  a_{22}& =b_{22}+c_{22}
\end{align}
$$


[^1]: Learn more: [LaTeX Equations](https://liascript.github.io/course/?https://raw.githubusercontent.com/liaScript/docs/master/README.md#52)


## Animations

Arrange content into sub-slides

    {{0-1}}
************

<article class="main-info">
Motor Trend Car Road Tests

The data was extracted from the 1974 Motor Trend US magazine, and comprises fuel consumption and 10 aspects of automobile design and performance for 32 automobiles.

</article>

***************

    {{1}}
***************
``` R
knitr::kable(head(mtcars)[,c("mpg",	"cyl", "disp", "hp", "wt")])
```
@runR

***************


## Incremental animations of list elements

Lists can optionally be displayed incrementally[^1]:

{1}{First item}

{2}{Second item}

{3}{Third item}


[^1]: Learn more: [Incremental animations](https://liascript.github.io/course/?https://raw.githubusercontent.com/liaScript/docs/master/README.md#103)

## Display various media 

You can also use the following as content elements[^1]

          {{0-1}}
**********************

-   An image

    - Syntax: `![Portrait of a lady](https://upload.wikimedia.org/wikipedia/commons/thumb/c/c3/Leonardo_da_Vinci_%28attrib%29-_la_Belle_Ferroniere.jpg/723px-Leonardo_da_Vinci_%28attrib%29-_la_Belle_Ferroniere.jpg "La Belle Ferronnière, c. 1490–1498")`

-   An image gallery

![Portrait of a lady](https://upload.wikimedia.org/wikipedia/commons/thumb/c/c3/Leonardo_da_Vinci_%28attrib%29-_la_Belle_Ferroniere.jpg/723px-Leonardo_da_Vinci_%28attrib%29-_la_Belle_Ferroniere.jpg "La Belle Ferronnière, c. 1490–1498")
![Lady with an Ermine](https://upload.wikimedia.org/wikipedia/commons/thumb/f/f9/Lady_with_an_Ermine_-_Leonardo_da_Vinci_-_Google_Art_Project.jpg/761px-Lady_with_an_Ermine_-_Leonardo_da_Vinci_-_Google_Art_Project.jpg "Lady with an Ermine, c. 1489–1491, Czartoryski Museum, Kraków, Poland")
![Mona Lisa](https://upload.wikimedia.org/wikipedia/commons/thumb/e/ec/Mona_Lisa%2C_by_Leonardo_da_Vinci%2C_from_C2RMF_retouched.jpg/687px-Mona_Lisa%2C_by_Leonardo_da_Vinci%2C_from_C2RMF_retouched.jpg "Mona Lisa or La Gioconda c. 1503–1516, Louvre, Paris")

*******************

          {{1}}
**********************

-   A video

    - YouTube: `!?[The Future of Programming](https://www.youtube.com/watch?v=8pTEmbeENF4)`

-   A sound clip 

    - Syntax: `?[soundcloud](https://soundcloud.com/glennmorrison/beethoven-moonlight-sonata)`

**************

[^1]: Learn more: [Media References](https://liascript.github.io/course/?https://raw.githubusercontent.com/liaScript/docs/master/README.md#17)


## Absolute Position

Position images or other elements at precise locations


![](https://quarto.org/docs/presentations/revealjs/demo/mini/images/kitten-400-350.jpeg)<!-- style="position: absolute; top: 170px; left: 30px; width: 400px; height: 400px;"-->

![](https://quarto.org/docs/presentations/revealjs/demo/mini/images/kitten-400-350.jpeg)<!-- style="position: absolute; top: 150px; right: 80px; width: 450px;"-->

![](https://quarto.org/docs/presentations/revealjs/demo/mini/images/kitten-400-350.jpeg)<!-- style="position: absolute; bottom: 110px; right: 130px; width: 300px;"-->


## Switch between graph and data table

    {{0-1}}
***************
``` R +Plot
png(file="out.png")
library(ggplot2)
ggplot(mtcars, aes(hp, mpg, color = am)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "loess")
```
@runR

**************

    {{1}}
***************

``` R +Data
knitr::kable(mtcars)
```
@runR

*************

## Quiz

Reading loud works best in Google Chrome and Edge. 

                    {{English Female |>}}
The film that I saw [[(that)|those|these|then]] night wasn’t very good.
It was all [[ about ]] a man [[ who ]] built a
time machine so he [[ could ]] travel back in time.
It took him ages and ages [[ to ]] build the machine.

[^1]: Learn more: [Quiz](https://liascript.github.io/course/?https://raw.githubusercontent.com/liaScript/docs/master/README.md#61)

## Display your organisation's logo
 
Use the icon macro to display your organisation's logo on top of the slides.

It works the same as for the logo macro.[^1]

```
<!--
# define this in the sectio on top of the markdown file
icon: ./pics/logo.png
-->
```

[^1]: Learn more: [Special macros](https://liascript.github.io/course/?https://raw.githubusercontent.com/liaScript/docs/master/README.md#176)


## Interactive Slides - widgets - show other templates

Include Jupyter widgets and htmlwidgets in your presentations

leaflet is not installed in the CodeRunner, most probably possible via templates

``` R
#| echo: false
#| fig-height: 5
library(leaflet)
leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addMarkers(lng=174.768, lat=-36.852, popup="The birthplace of R")
```
@runR

::: footer
Learn more: [Jupyter widgets](https://quarto.org/docs/interactive/widgets/jupyter.html), [htmlwidgets](https://quarto.org/docs/interactive/widgets/htmlwidgets.html)
:::

## Interactive Slides (probably out of scope / too R specific)

Turn presentations into applications with Observable and Shiny. Use component layout to position inputs and outputs.

```{r}
ojs_define(actors = data.frame(
  x = rnorm(100),
  y = rnorm(100)
))
```

```{ojs}
//| panel: sidebar
viewof talentWeight = Inputs.range([-2, 2], { value: 0.7, step: 0.01, label: "talent weight" })
viewof looksWeight = Inputs.range([-2, 2], { value: 0.7, step: 0.01, label: "looks weight" })
viewof minimum = Inputs.range([-2, 2], { value: 1, step: 0.01, label: "min fame" })
```

```{ojs}
//| panel: fill
import { plotActors } from './actors.js';
plotActors(actors, talentWeight, looksWeight, minimum)
```

::: footer
Learn more: [Observable](https://quarto.org/docs/interactive/ojs/), [Shiny](https://quarto.org/docs/interactive/shiny/), [Component Layout](https://quarto.org/docs/interactive/layout.html)
:::

## Easy Navigation and different views

Quickly jump to other parts of your presentation by opening the title panel on the left side. 

Presentation view, Slides view, and the Textbook view.

## Authoring Tools

- Live side-by-side preview for any text document in Github.dev

Learn more: 

- Or LiveEditor 

Learn more: [Jupyter](https://quarto.org/docs/tools/jupyter-lab.html), [VS Code](https://quarto.org/docs/tools/vscode.html), [Text Editors](https://quarto.org/docs/tools/text-editors.html)

## And More...

-   [Lia]() Android app (presentations look great on mobile, swipe to navigate slides)
-   [Footer & Logo](https://quarto.org/docs/presentations/revealjs/#footer-logo) (optionally specify custom footer per-slide)
-   [Auto-Slide](https://quarto.org/docs/presentations/revealjs/presenting.html#auto-slide) (step through slides automatically, without any user input)
-   [Multiplex](https://quarto.org/docs/presentations/revealjs/presenting.html#multiplex) (allows your audience to follow the slides of the presentation you are controlling on their own phone, tablet or laptop).

