<!--
title: "Liascript Presentations"

import: https://raw.githubusercontent.com/LiaScript/CodeRunner/master/README.md

link: style.css

@runR: @LIA.eval(`["main.R"]`, `none`, `Rscript main.R`)
-->

# Hello, There

This presentation will show you examples of what you can do with [Liascript](https://liascript.github.io), including:

-   Presenting code and LaTeX equations
-   Including computations in slide output
-   Image, video, and iframe backgrounds
-   Fancy transitions and animations
-   Printing to PDF

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

<article class="main-info">
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
</article>

<!-- class="sub-info" -->
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


## Incremental Lists

Lists can optionally be displayed incrementally:

-   First item
-   Second item
-   Third item


::: footer
Learn more: [Incremental Lists](https://quarto.org/docs/presentations/revealjs/#incremental-lists)
:::

## Auto-Animate {auto-animate="true" auto-animate-easing="ease-in-out"}

Automatically animate matching elements across slides with Auto-Animate.

::: r-hstack
::: {data-id="box1" auto-animate-delay="0" style="background: #2780e3; width: 200px; height: 150px; margin: 10px;"}
:::

::: {data-id="box2" auto-animate-delay="0.1" style="background: #3fb618; width: 200px; height: 150px; margin: 10px;"}
:::

::: {data-id="box3" auto-animate-delay="0.2" style="background: #e83e8c; width: 200px; height: 150px; margin: 10px;"}
:::
:::

::: footer
Learn more: [Auto-Animate](https://quarto.org/docs/presentations/revealjs/advanced.html#auto-animate)
:::


## Slide Backgrounds (check slide from online symposium)

Set the `background` attribute on a slide to change the background color (all CSS color formats are supported).

Different background transitions are available via the `background-transition` option.

::: footer
Learn more: [Slide Backgrounds](https://quarto.org/docs/presentations/revealjs/#color-backgrounds)
:::

## Media Backgrounds {background="#43464B" background-image="images/milky-way.jpeg"}

You can also use the following as a slide background:

-   An image: `background-image`

-   A video: `background-video`

-   An iframe: `background-iframe`

::: footer
Learn more: [Media Backgrounds](https://quarto.org/docs/presentations/revealjs/#image-backgrounds)
:::


## Auto-Animate {auto-animate="true" auto-animate-easing="ease-in-out"}

Automatically animate matching elements across slides with Auto-Animate.

::: r-stack
::: {data-id="box1" style="background: #2780e3; width: 350px; height: 350px; border-radius: 200px;"}
:::

::: {data-id="box2" style="background: #3fb618; width: 250px; height: 250px; border-radius: 200px;"}
:::

::: {data-id="box3" style="background: #e83e8c; width: 150px; height: 150px; border-radius: 200px;"}
:::
:::

::: footer
Learn more: [Auto-Animate](https://quarto.org/docs/presentations/revealjs/advanced.html#auto-animate)
:::

## Absolute Position

Position images or other elements at precise locations

    
![](https://quarto.org/docs/presentations/revealjs/demo/mini/images/kitten-400-350.jpeg)<!-- style="top: 170px; left: 30px; width: 400px; height: 400px;"-->

![](mini/images/kitten-450-250.jpeg){.absolute .fragment top="150" right="80" width="450"}

![](mini/images/kitten-300-200.jpeg){.absolute .fragment bottom="110" right="130" width="300"}

::: footer
Learn more: [Absolute Position](https://quarto.org/docs/presentations/revealjs/advanced.html#absolute-position)
:::

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


## Interactive Slides - widgets - show other templates

Include Jupyter widgets and htmlwidgets in your presentations

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

## Interactive Slides {.smaller transition="slide"}

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

## Preview Links

Navigate to hyperlinks without disrupting the flow of your presentation.

Use the `preview-links` option to open links in an iframe on top of your slides. Try clicking the link below for a demonstration:

::: {style="text-align: center; margin-top: 1em"}
[Matplotlib: Visualization with Python](https://matplotlib.org/){preview-link="true" style="text-align: center"}
:::

::: footer
Learn more: [Preview Links](https://quarto.org/docs/presentations/revealjs/presenting.html#preview-links)
:::

## Adapt your icon

10 Built-in Themes (or [create your own](https://quarto.org/docs/presentations/revealjs/themes.html#creating-themes))

::: {layout-ncol="2"}
![](images/moon.png)


## Easy Navigation

::: {style="margin-bottom: 0.9em;"}
Quickly jump to other parts of your presentation
:::


::: footer
Learn more: [Navigation](https://quarto.org/docs/presentations/revealjs/presenting.html#navigation-menu)
:::

## Presentation view

::: {style="margin-bottom: 0.9em;"}
Free form drawing and slide annotations
:::



## Textbook view

Press `o` to toggle overview mode:

![](images/overview-mode.png){.border}

Hold down the `Alt` key (or `Ctrl` in Linux) and click on any element to zoom towards it---try it now on this slide.

::: footer
Learn more: [Overview Mode](https://quarto.org/docs/presentations/revealjs/presenting.html#overview-mode), [Slide Zoom](https://quarto.org/docs/presentations/revealjs/presenting.html#slide-zoom)
:::

## Slides view

Press `s` (or use the presentation menu) to open speaker view

![](images/speaker-view.png){fig-align="center" style="border: 3px solid #dee2e6;" width="780"}

::: footer
Learn more: [Speaker View](https://quarto.org/docs/presentations/revealjs/presenting.html#speaker-view)
:::

## Authoring Tools

Live side-by-side preview for any text document in Github.dev

::: columns
::: {.column width="50%"}
![](images/jupyter-edit.png){.border .border-thick}
:::

::: {.column width="50%"}
![](images/jupyter-preview.png){.border .border-thick}
:::
:::

::: footer
Learn more: [Jupyter](https://quarto.org/docs/tools/jupyter-lab.html), [VS Code](https://quarto.org/docs/tools/vscode.html), [Text Editors](https://quarto.org/docs/tools/text-editors.html)
:::

## Authoring Tools

LiveEdit mode

![](images/rstudio.png){.border width="900"}

::: footer
Learn more: [RStudio](https://quarto.org/docs/tools/rstudio.html)
:::

## And More...

-   [Touch](https://quarto.org/docs/presentations/revealjs/advanced.html#touch-navigation) optimized (presentations look great on mobile, swipe to navigate slides)
-   [Footer & Logo](https://quarto.org/docs/presentations/revealjs/#footer-logo) (optionally specify custom footer per-slide)
-   [Auto-Slide](https://quarto.org/docs/presentations/revealjs/presenting.html#auto-slide) (step through slides automatically, without any user input)
-   [Multiplex](https://quarto.org/docs/presentations/revealjs/presenting.html#multiplex) (allows your audience to follow the slides of the presentation you are controlling on their own phone, tablet or laptop).

