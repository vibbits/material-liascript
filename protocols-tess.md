<!--
author:   Alexander Botzki

email:    alexander.botzki@vib.be

version:  0.1

language: en

icon: img/elixir-tess-b37789ea9aced77e795a3c2dc3a2fa583d9dc0e0eb392d589a8cd8d534be8067.svg

narrator: US English Female

comment:  the course is based on the publication.

link:     https://cdn.jsdelivr.net/chartist.js/latest/chartist.min.css

script:   https://cdn.jsdelivr.net/chartist.js/latest/chartist.min.js

link: style.css

@orcid: [@0](@1)<!--class="orcid-logo-for-author-list"-->

-->

# Protocols for efficient use of TeSS

Hello and welcome, this is our interactive hands-on course about efficient use of the ELIXIR TeSS platform.

> We are using the interactive Open Educational Resource online/offline course infrastructure called LiaScript. 
> It is a distributed way of creating and sharing educational content hosted on github.
> To see this document as an interactive LiaScript rendered version, click on the
> following link/badge:
>
> [![LiaScript](https://raw.githubusercontent.com/LiaScript/LiaScript/master/badges/course.svg)](https://liascript.github.io/course/?https://raw.githubusercontent.com/vibbits/material-liascript/master/protocols-tess.md)

To prepare this training material, we have reused the content of the publication 'Making Bioinformatics Training Events and Material More Discoverable Using TeSS, the ELIXIR Training Portal' published in Current Protocols in Bioinformatics. The original publication of the protocols can be found [here](https://doi.org/10.1002/cpz1.682). 

--------------------------------------------

*About ELIXIR Training Platform*

The ELIXIR Training Platform was established to develop a training community that spans all ELIXIR member states (see the list of Training Coordinators). It aims to strengthen national training programmes, grow bioinformatics training capacity and competence across Europe, and empower researchers to use ELIXIR's services and tools. 

One service offered by the Training Platform is TeSS, the training registry for the ELIXIR community. Together with ELIXIR France and ELIXIR Slovenia, VIB as lead node for ELIXIR Belgium is engaged in consolidating quality and impact of the TeSS training resources (2022-23) (https://elixir-europe.org/internal-projects/commissioned-services/2022-trp3).

The Training eSupport System was developed to help trainees, trainers and their institutions to have a one-stop shop where they can share and find information about training and events, including training material. This way we can create a catalogue that can be shared within the community. How it works is what we are going to find out in this course.

*About VIB and VIB Technologies*

VIB is an entrepreneurial non-profit research institute, with a clear focus on groundbreaking strategic basic research in life sciences and operates in close partnership with the five universities in Flanders – Ghent University, KU Leuven, University of Antwerp, Vrije Universiteit Brussel and Hasselt University.

As part of the VIB Technologies, the 12 VIB Core Facilities, provide support in a wide array of research fields and housing specialized scientific equipment for each discipline. Science and technology go hand in hand. New technologies advance science and often accelerate breakthroughs in scientific research. VIB has a visionary approach to science and technology, founded on its ability to identify and foster new innovations in life sciences.

The goal of VIB Technology Training is to up-skill life scientists to excel in the domains of VIB Technologies, Bioinformatics & AI, Software Development, and Research Data Management.

--------------------------------------------

*Editorial team for this course*

Authors: @[orcid(Alexander Botzki)](https://orcid.org/0000-0001-6691-4233), @[orcid(Bruna Piereck)](https://orcid.org/0000-0001-5958-0669)

Contributors: Finn Bacall, Aitor Apaolaza, Munazah Andrabi, Chris Child, Carole Goble, Olivier Sand

Technical Editors: Alexander Botzki

License: [![CC BY](img/picture003.jpg)](http://creativecommons.org/licenses/by/4.0/)

## Introduction

                {{0}}
*************************

Many trainers and their organisations are passionate about sharing their training material. Sharing training material has several benefits like providing a record of recognition as an author, offering inspiration to other trainers, enabling researchers to discover training resources for their personal learning path, and facilitating to complement the training resource landscape by data-driven gap analysis from the bioinformatics community (see [PLOS Computational Biology, 16(5), e1007854](https://doi.org/10.1371/journal.pcbi.1007854)).

![advantages of TeSS](img/XgXcCp3.png)

In this tutorial, we elaborate a series of protocols leading the TeSS users through the process of searching and filtering TeSS content (Scenario 1 to 3). Scenario 4 to 7 show how to register manually and automatically training events and material in the [Training eSupport System (TeSS), ELIXIR’s training portal](https://tess.elixir-europe.org). More information about TeSS can be found in [this publication](https://doi.org/10.1093/bioinformatics/btaa047). 

TeSS provides a one-stop shop for trainers and trainees to discover online information and content, including training materials, events and interactive tutorials. Following these protocols, you will contribute to promote your training events and add to a growing catalog of materials. This will concomitantly increase the FAIRness of training materials and events.

Training registries like TeSS use a scraping mechanism to aggregate training resources from many providers when they have been annotated using training specific [BioSchemas specifications](https://doi.org/10.1101/2022.11.24.516513). Protocol 8 will teach you how to enrich your training resources to allow for more efficient sharing of the structured metadata like prerequisites, target audience, and learning outcomes using Bioschemas specification. Once more and more training events and material are aggregated in TeSS, searching the registry for specific events and material becomes crucial.

![overview annotations](img/KbCIsol.jpg)

In the next section, you will find links where you can follow all scenarios.

**************************

                     {{1}}
***********

[Scenario 1](#4): Searching for training events in TeSS

[Scenario 2](#5): Integrating TeSS widgets on your web site

[Scenario 3](#6): Searching for e-learning materials in TeSS

[Scenario 4](#7): Log into TeSS using an institutional account

[Scenario 5](#8): Manual registration of a training event in TeSS

[Scenario 6](#9): Manual registration of training material in TeSS

[Scenario 7](#10): Registration of a content provider in TeSS

[Scenario 8](#11): Automated harvesting of a training event in TeSS

***********

## More background information about TeSS

During the ELIXIR-EXCELERATE project, the TeSS platform has been established as the key application for bioinformatics training events and material aggregation. During a period of further adoption of the TeSS platform as such on the international level e.g. in Australia and in the Netherlands, it became obvious that the implementation of individual automatic content scrapers for more than 30 content providers lead to an enormous diversity of the code base and finally to a maintenance challenge. Therefore, during the last two to three years, a strategic reorientation of the automated harvesting procedure has been proposed and implemented in TeSS using the Bioschemas profiles for Course, CourseInstance and TrainingMaterial. This approach will put the content providers in a more prominent position since they are in turn responsible for a more structured delivery of more comprehensive metadata of their training events and material to TeSS as the aggregation platform. In the longer run, this will increase the quality level of the aggregated resources and its usefulness to TeSS user base.
TODO: Would be insteresting to add a box explaining the specificity of each?

Content providers looking to promote use of particular tools and databases via their course offering or the specific training material can find information about them in other ELIXIR registries, such as [bio.tools](https://bio.tools) for tools and web services  and [FAIRsharing.org](https://fairsharing.org) for databases, standards and policies . Lately, such registries started to track training relating to the resources they list. Therefore, it is mutually beneficial that TeSS provides the relevant links, allowing users to search for resource-specific training; and the resource pages in bio.tools and FAIRsharing that contain reciprocal links for users to discover relevant training in TeSS.

TODO: Are TeSS, BioTools and FAIRsharing connected? In content? 

Related to the provisioning of persistent identifiers, TeSS provides the Digital Object Identifier (DOI) as a unique means to identify training material. Given that minting services by e.g. university libraries are not yet broadly available for teaching and learning material in contrast to research datasets and collections, associated workflows, software and models, we recommend to use services like OSF or zenodo to get own unique, persistent URLs and DOIs.

Given the engagement of the ELIXIR Train the Trainer instructor's community to create a vibrant environment for reciprocal support and exchange of experiences via the ELIXIR Train the Trainer programme, the TeSS content providers are encouraged to  follow best-practice guidance on course and training material development (Via et al., 2020). Therefore, we recommend applying the Bloom hierarchy of cognitive skills to formulate the Learning objectives for the courses as well as for training material. For more detailed information about our pedagogical model, please browse to the course material of the [ELIXIR Train the Trainer course](https://github.com/TrainTheTrainer/ELIXIR-EXCELERATE-TtT). 


## Scenario 1: Searching for training events in TeSS

            {{0}}
*********************

Training resources (both events and materials) may be searched in TeSS in several ways. If you are on the main page of TeSS, a general search can be performed based on keywords, which will return separate lists of events and materials. Alternatively, events or materials can be searched for independently of each other. This second approach allows more precise filtering on several parameters (e.g., event type, country, and target audience) alone or in combination.

![screenshot for TeSS page with search box](img/GxGphHo.png)

In this protocol, we provide examples of searches for specific events or materials. Be aware that by searching for events, the standard search returns only future events. If you want to include past events, you need to click the **Show past event** button. 

**********************

              {{1}}
**********************

Please navigate to [TeSS](https://tess.elixir-europe.org) and follow the steps shown below.

1. Select Events in the top menu.
2. **Search for** ‘Single-cell’ as the keyword in the search text box on the right side of the page. **You can click on the magnifier icon or press Enter**.
3. Click on the **Show past events** button **at the bottom of the side menu**. Click on 'Filters' if you do not see this side menu due to your browser being non-fullscreen.
4. On the side menu, scroll down to the Country facet. We will select the United Kingdom.
5. On the side menu, scroll to the facet Target audience. Select ‘Graduate students’.
6. As a result, you get the list of training events about single-cell data analysis available in the United Kingdom for Graduate students.
7. Change the sort filter to 'Latest' in order to see an anti-chronological listing of the events.

At any time you may remove any of the filters you applied to reduce stringency of your search.

Please proceed to the next step to see the results.

*********

             {{2}}
**********************

At any time you may remove any of the filters you applied to reduce stringency of your search.

In June 2023, we found 13 past events with this search query.  

Have a look at the URL to understand how you could query TeSS via the API.

https://tess.elixir-europe.org/events?country=United+Kingdom&include_expired=true&q=single-cell&target_audience=Graduate+students&sort=late

Let's have a look at the documentation of the API: https://tess.elixir-europe.org/api/json_api#tag/events

Here you see all the query parameters. In our query, we have been using several of these parameters and concatenated them.

1. First, we used the free text query parameter `q` and set it to `single-cell`. 
2. Second, we used the `include_expired` parameter and set it to `true`. 
3. Third, we used the `country[]` parameter and set it to "United+Kingdom". Note that you need to replace the space by a `+` sign. 
4. Fourth, we used the `target_audience[]` parameter and set it to `Graduate+students". 
5. To change the sorting, use the `sort` parameter. You have 5 options, we have used `late` to sort by date, latest to earliest (events only).  

If you would like to get the events from France, too, you'd need to add yet another query with `country[]=France`.

To get the full query, you need to concatenate all the individual queries with an ambersand `&`. 

*********

## Scenario 2: Integrating TeSS widgets on your web site

            {{0}}
*********************

Widgets are chunks of JavaScript code that can be copied into website source code to display TeSS content. There are several different styles and functionalities available from our configurable widget suite. Widgets can be used to enhance your site and offer your community lists of relevant events or training resources. There are examples of code on [this page](https://elixirtess.github.io/TeSS_widgets/).

********************

            {{1}}
*********************

You can use these widgets to enhance your site and offer your community lists of relevant events or training resources! Check the examples of code on the [widget examples page](https://elixirtess.github.io/TeSS_widgets/). We will provide one example below.

**How to implement a tabular layout:**

1. Choose the specific widget type **Events table with a sidebar**. A working example is shown on the github pages from the ELIXIR TeSS team at https://elixirtess.github.io/TeSS_widgets. 

2.	Manually annotate the HTML web site e.g. index.html with the code snippet for the TeSS widget. More specifically, enter the code snippet as shown below as a sibling of a div element.

```html
<link rel="stylesheet" property="stylesheet" href="https://elixirtess.github.io/TeSS_widgets/css/tess-widget.css"/>
<div id="tess-widget-events-table" class="tess-widget tess-widget-faceted-table"></div>
<script>
function initTeSSWidgets() {
    TessWidget.Events(document.getElementById('tess-widget-events-table'),
        'FacetedTable',
        {
            opts: {
                columns: [{name: 'Date', field: 'start'},
                    {name: 'Name', field: 'title'},
                    {name: 'Location', field: 'location'}],
                allowedFacets: ['scientific-topics', 'country', 'city', 'target-audience'],
                facetOptionLimit: 5
            },
            params: {
                q: 'Python',
                country: ['Belgium', 'United Kingdom']
            }
        });
}
</script>
<script async="" defer="" src="https://elixirtess.github.io/TeSS_widgets/js/tess-widget-standalone.js" onload="initTeSSWidgets()"></script> 
```

3.	Save the HTML page and serve the HTML page for visual inspection.

4.	Validate visually that the rendered HTML page includes a table with the Events extracted from TeSS. It will show a filtered list according to the values in the params section (see snippet from step 2). In this case, all future courses aggregated in TeSS will be live filtered on the general search term ‘Python’ and display only the courses from ‘Belgium’ and ‘United Kingdom’.

5. You can adapt the code:
- Remove the term ‘Python’ from the general search. Now, the table underneath should be updated with all the up-coming courses from Belgium and the United Kingdom. 
- By adapting other facets like ‘Scientific Topics’ and ‘Target audience’, other filtering schemes can be established. 

At any time you may remove any of the filters you applied to reduce stringency of your search. If you would like to include past events in your search, click on the Show past events button.

If you'd like to use these widgets for training materials, investigate the example for training materials on the [same website](https://elixirtess.github.io/TeSS_widgets/).

In the next step, we will show how to VIB Technology Training has implemented the widgets on their home page.


*****************

            {{2}}
*********************

TODO: find an events widget

If you'd like to see an implementation of the TeSS widget for materials, please browse to [https://technologytraining.sites.vib.be/en/overview-offering](https://technologytraining.sites.vib.be/en/overview-offering) and scroll to the bottom of the page.

*********************

## Scenario 3: Searching for e-learning materials in TeSS

            {{0}}
*********************

E-learning resources are a specific type of training materials in TeSS which are displayed on the dedicated section **E-learning** on the main TeSS page. 

![screenshot for TeSS page with search box](img/CXt2MEE.jpg)

But you can also perform a specific search on the Materials page filtering on 'Resource type' and retrieve e-learning material by specifying the value 'e-learning'.

*********************

            {{1}}
*********************

Please navigate to [TeSS](https://tess.elixir-europe.org).

1.	Select Materials in the top menu.

2.	Enter a keyword in the search text box on the right side of the page. We will enter ‘Data retrieval’.

3.	On the side menu, select a difficulty level. We will select ‘Beginner’.

4.	Scroll to the facet ‘Resource type’ and select ‘e-learning’.

5.	You get the list of training materials about data retrieval intended for a beginner level audience whereby the resource has been tagged as e-learning material. At any time you may remove any of the filters you applied to reduce stringency of your search.

********************

## Scenario 4: Log into TeSS using an institutional account

            {{0}}
*********************

In order to be able to register resources in TeSS you need to log in the registry. The following steps are specific for the login procedure with the Life Sciences Account. Alternatively, you can use another existing account like Google, Apple or an ORCID. In this case you will be redirected to the respective login portal.

**We will not be able to provide a detail procedure for all the federated authentication mechanisms here. Try the example below.**

*********************

            {{1}}
*********************

1.	Open your browser and go to [TeSS](https://tess.elixir-europe.org/). 

2.	Click the “Log In” dropdown button on the top-right corner of the page.

3.	To log in using your institutional account, choose Log in with LS Login from the dropdown menu
You will be redirected to the LifeScience RI authentication page (Fig. 1).

![screenshot for login page](img/6dXdfgu.jpg "Figure 1")

4.	Start typing the name of your own institution, for instance ‘VIB’, into the text box then choose the appropriate option that appears and proceed with your usual institutional login procedure (Fig. 2). 

![screenshot for login page](img/l5mQUND.png "Figure 2")

5. You should now be logged into TeSS. Once successful, you should be taken to the TeSS Welcome page with a message stating “Logged in successfully.”  You should also see your username in the top bar of the page, which you can click to view and edit your TeSS profile.

![screenshot for login page](img/meCcJjO.png "Figure 3")

6.	Check all the fields in your user profile are correct and click the “Update Profile” button.
You will also see a button to log out of TeSS.

*********************

## Scenario 5: Manual registration of a training event in TeSS

            {{0}}
*********************

TODO: move to Scenario 7?
A training event in TeSS is described by the BioSchemas profile CourseInstance, while the schema.org Event profile is something different. So, let’s start this session pointing out a few differences of each profile:

An unique event such as a congress, symposium, webinar will be classified using the Event profile in Schemas.org, while the Bioschemas Course profile is attached to a unit of teaching that typically lasts one or more academic terms, is led by one or more instructors (teachers or professors), and has a fixed roster of students. A course usually covers an individual subject. Such a course will be linked to several course instances that happen at different moments, locations and might be also subjected to updates and small differences in the material. A course and/or course instance should be linked to a course material, that represents the content of this course/course instance. In TeSS, though, courses and courses instances are in called Event in the interface, while this is not the same as the profile in schemas.org. It works to filter among all courses, webinars, symposia, congresses etc.

Such an event in TeSS is a link to a single training event sourced by a provider along with description and related meta information (e.g. date, location, audience, ontological categorization, keywords, etc.). 

Training events in TeSS can be added manually or automatically harvested from a provider's website. Here is the procedure to register for training events manually.

![screenshot for login page](img/rhh0rv5.png)

*********************

            {{1}}
*********************

1. Login to TeSS at https://tess.elixir-europe.org following [Scenario 3](#6).

2. To register for a training event, click on the Events menu.

3. Click on “+ Register event”. As an example, we use an RNAseq training event.
You will be directed to a form type of page. 

![screenshot for login page](img/appDpje.png)

4. Select that Type at the beginning of the form since we enter a face-to-face event.
If your event is online, choosing the option "Online" which will modify the form to fit that type of event i.e. no address field has to be filled in. 

5. In the Title field, enter the Title of the event, which is “Bulk RNASeq: from counts to differential expression” in our example.

6. The optional Subtitle field can be filled in case you’d like to specify more context on the event, like “Spring edition” of a regularly repeated course.

7. In the URL field, enter the URL where the announcement and registration of the training event can be found. In our example: “https://training.vib.be/all-trainings/bulk-rnaseq-counts-differential-expression-2”

8. In the subsequent field Description, provide general information and relevant context about the training event in the form of short text, like in our example: “The course consists of introductory online material on producing count matrices and two face-to-face sessions on differential expression analysis in R and all the questions that arise when trying the analysis on your own data.”
This field supports markdown syntax.

********************

            {{2}}
*********************


9. Specify the start date and door time of the event in the field Start via the pop-up Date and Time Picker. In our example, the event starts on the 29th of September 2022.

10. Specify the end date of the event in the field End via the pop-up Date and Time Picker. Our example event will end on the 4th of October 2022.

11. In the field Timezone, enter the time zone of the location where the event will occur and pick the one that applies. In our example, we select ‘(GMT:+01:00) Amsterdam’ since the event takes place in Leuven, Belgium.

12. In the field Duration, indicate how long the event will last in hours, days, whatever applies best. Enter ‘2 days’ for our example.

13. In the field Prerequisites, enter which prerequisites are necessary to follow the training event. For example, write “Background knowledge of NGS data formats and the first steps in the analysis workflow (fastqc -> bam files). If you are a newbie in the field, you have to follow the NGS introduction training first. Experience in R programming is highly recommended. If you never worked with R, you should attend the R introduction training first.”
This field supports markdown syntax.

14. In the field Learning objectives, list the learning objectives of the training event like 
○	Learn how RNA-seq reads are converted into counts
○	Understand QC steps that can be performed on RNA-seq reads
○	Generate interactive reports to summarise QC information with MultiQC
○	Use the Galaxy Rule-based Uploader to import FASTQs from URLs
○	Create a Galaxy Workflow that converts RNA-seq reads into counts
We recommend applying the Bloom hierarchy of cognitive skills to formulate the Learning objectives. For more detailed information about our pedagogical model, please browse to the course material of the ELIXIR Train the Trainer course https://github.com/TrainTheTrainer/ELIXIR-EXCELERATE-TtT This field supports markdown syntax. 

15. For the Address and related fields, enter the relevant information as
○	Address: fill in the address by starting to type the address. Once you select an address, the Google maps view will be shown and the fields Venue, City, Region, Country and Postcode will be filled in automatically
○	Venue: the venue of the event e.g. street
○	City: the city where the venue is located
○	Region: the region where the city is located
○	Country: the country where the event will take place
○	Postcode: the postcode of the venue
In our example, these are
○	Herestraat 49
○	Campus Gasthuisberg
○	Leuven
○	Flemish region
○	Belgium
○	3000

*******************

            {{3}}
*********************


16. For the Eligibility field, there are three options you can choose from: First come first served, Registration of interest, By invitation. You could provide more than one eligibility type, if applicable. In our example, we enter First come, first serve.

17. Enter for example ELIXIR-Belgium in the Organiser field.
This field is a free text field.

18. Enter the name of the training coordinator in the contact field. 
Since this is also a free text field, one or more contact persons preferably with an email address can be added.

19. Start typing ‘VIB’ as the name of the host institution(s) in the ‘Host institutions’ field. 
You will get a dropdown menu to choose from. If your institution is not present in the menu, add it. Add as many as you’d like.

20. Enter ‘transcriptomics’ in the Keywords field.
You can start typing keyword(s) and you will get a list to choose from. Add as many as you want.

21. Enter ‘Life scientists researchers’ as Target audience. 
Start typing target audience types and you will get a list to choose from. Add as many as you want.

22. In the field Scientific topics, enter ‘RNA-Seq, Sequencing’.
You can enter the appropriate scientific topics according to the EDAM ontology. In case you do not have a matching entry, you could use the keywords field.

23. Enter ‘Differential gene expression profiling’ in the Operations field.
you can enter the appropriate operations according to the EDAM ontology. In case you do not have a matching entry, you could use the keywords field.

*******************

            {{4}}
*********************

24. Add e.g. 20 as maximum number of attendants to the event in the field Capacity. 

25. Enter ‘Workshops and courses’ as Event type.
There are four options for the event type: Workshops and courses, Awards and prizegivings, Meetings and conferences, Reception and networking. You could provide more than one event type if applicable.

25. As Tech requirements, enter “R, Rstudio, various R packages”. 
You can list here any technical requirement needed to follow the training event. Ideally, you also provide a link to the installation instructions. This field supports markdown syntax.

26. Enter ‘Certificate of attendance recognised by Doctoral School’ as the type of credit or recognition of attendance that will be delivered.

27. In the field External resources, search for ‘DESeq2’ in the tools section. Once the tool has been found in bio.tools, click on the “+” sign to add it to the entry.
You can also associate policies, standards and databases (provided by FAIRsharing) with your training event. In each case, start typing keywords and a list of existing resources will be listed to select from (click on the “+” sign).

28. Select ‘Cost to non-members’ as Cost basis.
There are three choices: Free to all, Cost to non-members, Cost incurred by all.

29. Select ‘VIB Training’ as Content Provider. 
In case you do not find your organisation, please register it first as a Provider in TeSS. 

**************

            {{5}}
*********************

30. Search for ‘transcriptomics data analysis made easy’ in the Materials search box and select suggested training materials to associate with this course. Start typing keywords and a list of existing training materials will be listed to select from (click on the “+” sign).

31. Select ‘ELIXIR Belgium’ as a node from the downtown list. Select the node the organiser(s) belong(s) to.

32. Enter ‘VIB’ as Sponsors/Funders of your training event.

33. Click on the Add event button to finalize the manual registration.

34. The resulting course can be visited [here](TODO).

*******************

## Scenario 6: Manual registration of training material in TeSS

In the context of TeSS, a training material is a link to a single online training material sourced by a content provider (such as a text on a Web page, presentation, video, etc.) along with description and related meta information (e.g. ontological categorization, keywords, etc.). Materials can be added manually or automatically harvested from a provider's website. Here is the procedure to register training materials manually.

![sketch of components of training material](img/MOGjyQj.png)

1.	As outlined in [Protocol 4](#4), in order to register for training material in TeSS, login first. 

2.	To register for training material, click on the Material menu.

![screenshot about registering training material in TeSS](img/uKQ8Xr9.jpg)

3.	Click on “+ Register training material”. You will be directed to a form. We will use training material about the use of AlphaFold for structure prediction of proteins and binary protein complexes.

4.	Enter your title of the training material e.g. ‘How to use AlphaFold for prediction of single proteins or protein complexes using the HPC’.

5.	Enter the URL where the training material including the description can be found. For the above example e.g. ‘https://elearning.bits.vib.be/courses/alphafold/’

6.	Enter a description of the course in the next step. 

For our example, it would be:
“Architectural details, code and trained AlphaFold models were released by DeepMind in 2021. Given the high computational cost of deep learning algorithms, specialized hardware and software are required. Online solutions are available but come with considerable disadvantages. Therefore, the Flemish Supercomputer Center (VSC) provides high performance computing facilities, on which AlphaFold is installed and fully operational. This course gives a solid introduction on how AlphaFold can be easily and swiftly accessed using the HPC.
This tutorial material was created by Jasper Zuallaert (VIB-UGent), with the help of Alexander Botzki (VIB) and Kenneth Hoste (UGent). For questions and remarks, feel free to contact jasper.zuallaert@vib-ugent.be“

This field supports markdown syntax
                  
7.	In the field Keywords, start typing keyword(s) and you will get a list to choose from. Add as many as you want. New ones can also be added.

For our example: ‘Structure prediction’, ‘AlphaFold Database’, ‘Protein complex prediction’.

8.	Select the license of the material.We advise to use ‘Creative Commons Attribution Non commercial Share Alike 4.0 International’.
From the dropdown list, choose the one that applies. The dropdown includes “License not specified” in case you do not hold one for your material.

9.	Select the status that fits the material you’d like to register. 
There are three stati available: Active, Under development, Archived.

10.	Enter your contact email address in the Contact field. Enter one of more contact persons preferably with an email address.

11.	Enter DOI if available. In our example, we do not provide a DOI.
Add the URL of the DOI in this field if you have created a DOI for your training material. If you host your material in github, you can get a DOI by linking the github repository to Zenodo.

12.	You can provide a version number for the material, if available. In our case, we specify ‘1.0’ as in the Version field. 

13.	Please select the content provider from the dropdown list. In our case, it is ‘VIB Training’ as Content Provider.

14.	Enter the name of authors of the material. Here we write ‘Jasper Zuallaert’ and ‘Alexander Botzki’ in the Authors field. 

15.	Enter the name of contributors to the material.Preferentially, specify your ORCID. In our example, the name of the colleague from the HPC team ‘Kenneth Hoste’ is added in the Contributors field.

16.	Enter the event which is registered in TeSS where the training material is used. We recommend doing the linking from the event side (when the event is registered) according to Basic Protocol 4 step 30. In our example, it is‘Structural Prediction of Proteins using AlphaFold on the HPC’.

17.	Add the Target Audience. Start typing target audience types and you will get a list to choose from. Add as many as you want. In our case, we have entered ‘Life scientists with programming skills’ as Target audience.

18.	If prerequisites are necessary to follow the training event, list them here. This field supports markdown syntax. 
We added the following text as Prerequisites in the field: ‘You are encouraged to use your own laptop. For those who do not have a laptop, the YASARA software can be run in a remote Linux environment (access to cloud via web browser). Knowledge of command line and basic Python skills are recommended.’ 

19.	You can choose between three levels: Beginner, Intermediate, Advanced for the competency level. Our course is on an ‘intermediate’ level. 

20.	We recommend applying the Bloom hierarchy of cognitive skills to formulate the Learning objectives. For more detailed information about our pedagogical model, please browse to the course material of the ELIXIR Train the Trainer course https://github.com/TrainTheTrainer/ELIXIR-EXCELERATE-TtT. This field supports markdown syntax. 
 
![overview of the Bloom taxonomy](img/tl5eM8b.png)

In our example, we enter the following text as Learning Objectives: ‘Understand the technical methodology of AlphaFold2, understand the technical setup at the Flemish SuperComputer, Predict three-dimensional protein models with AlphaFold2 using the HPC at the VSC UGhent’

21.	Pick the creation data of your material as 'Date created'.

22.	Enter no value in Date modified. 
Pick the date the material was last modified, if applicable.

23.	Enter the date when the course has been published. 

24.	Enter the media used for your material (e.g. video, slides, pdf…). Just start typing a media type and you will get a list to choose from.
In our case, the Resource types are ‘Tutorial’ and ‘Slides’ .

25.	You can enter the appropriate scientific topics according to the EDAM ontology. In case you do not have a matching entry, you could use the keywords field (see above). 

In our example, we associate ‘Structure prediction, Protein Structure Analysis, Machine Learning’ as Scientific topics. 

26.	Click on ‘Suggested tools to associate with this resource’ in the External resources field.
You can associate policies, standards and databases (provided by FAIRsharing) or tools (provided by bio.tools) with your training event. In each case, start typing relevant words and a list of existing resources will be listed to select from (click on the “+” sign).

In our example, we enter ‘AlphaFold’ in the search field. Wait a couple of seconds until the result from the request to the ELIXIR Tools and Services Registry is shown. Select the entry ‘AlphaFold 2’ by clicking on the “+” sign.

27.	Select the node(s) from which the organiser(s) belong(s) to from the dropdown list. Here is it ‘Belgium’ which is selected in the Nodes field. 

28.	In the next step, you can enter the appropriate scientific topics according to the EDAM ontology. In case you do not have a matching entry, you could use the keywords field (see above). 

For our example, we select ‘Multiple Sequence Alignment, Structure Visualisation’ in the field Operations. 

29.	Click on ‘Register training material’.

30. The resulting training material can be visited in TeSS. 

![screenshot for login page](img/9COCE0h.png)

Our example can be found [here](https://tess.elixir-europe.org/materials/prediction-of-protein-structures-and-complexes-with-alphafold-on-the-hpc).

## Protocol 7: Registration of a content provider in TeSS

            {{0}}
*********************

Training resources (events and materials) may be added to TeSS to reach bigger audiences, increase impact and bolster event attendance. Registering events and training materials makes them more findable in a variety of ways to various user bases. TeSS features content providers which are entities (such as academic institutions, non-profit organisations, portals etc.) that provide training materials of relevance to life sciences and ELIXIR. In order to have a training event automatically harvested in TeSS, as described in [Scenario 8](#9), a content provider has to be registered first. It also adds visibility to content providers. Here is the procedure to register a new content provider in TeSS.

*********************

            {{1}}
*********************

1. Login to TeSS at https://tess.elixir-europe.org following [Scenario 4](#4).

2. To register a content provider, click on the Providers menu. 

3. Click on “+ Register content provider”. You will be directed to a form type of page.

4.  As an example, we will create the content provider NanoCommons. Enter “NanoCommons“ in the Title field.

5. In the URL field, enter the URL of the content provider, which is “https://www.nanocommons.eu/” in our example.

6. Enter “Jean Dupont” as a contact person in the Contact field.

7. In the subsequent field Description, provide general information and relevant context about the content provider in the form of short text.

This field supports markdown syntax.

Here is one example: “NanoCommons will deliver a sustainable and openly accessible nanoinformatics framework (knowledgebase and integrated computational tools, supported by expert advice, data interpretation and training), for assessment of the risks of NMs, their products and their formulations. NanoCommons combines Joint Research Activities to implement the nanoinformatics Knowledge Commons, Networking Activities to facilitate engagement with the research community, industry and regulators, and provision of funded Access to the nanoinformatics tools via funded calls for Transnational Access.” 

***************

            {{2}}
*********************

8. Enter a image to have a **visual representation of the content provider**. The image can be added by URL of the image location or by uploading an image.

![screenshot for login page](img/qyOSDyb.png)

9. Enter the type of the content provider in the Type field. You can choose ‘Project’, ‘Organisation’ or ‘Portal’.

10. In the field Approved Editors, You can add more than one **registered** user as Approved Editor. Selected users can be removed by clicking on the red cross.

11. In the field Keywords, enter the keywords related to the content provider.

12. Select an associated ELIXIR node if applicable. In this example, there is no direct association with ELIXIR nodes, so we leave it empty.

13. Confirm the creation of the content provider by clicking on ‘Register content provider’.

******************

## Scenario 8: Automated harvesting of a training event in TeSS

            {{0}}
*********************

Registering events and training materials makes them more findable in a variety ways to various user bases. TeSS features several options for automatically "harvesting" resources from external sources. This can be helpful if you are maintaining a large collection of events and materials that changes frequently. 

![screenshot for login page](img/9COCE0h.png)

To register resources in TeSS automatically, we need to be able to extract data from target sources reliably. To this end, it is helpful if the data is structured according to a globally used standardized format. 

The following are examples of the kinds of structured data that TeSS can work with:

If your website currently includes no structured data, and you’d like your resources added to TeSS, we recommend using Bioschemas to structure your site. [Schema.org](https://schemas.org) is a project running by a consortium of search engines that has created an extensive library of schemas. Web-masters can use these schemas to explicitly mark-up their websites content in order to improve search engine visibility and interoperability. 

[Bioschemas](https://bioschemas.org) is an initiative to supplement the work of schema.org to help improve the findability of online resources in the life sciences. 

TeSS supports the following Bioschemas profiles: for events that are courses **CourseInstance** and **Course**, for other events **Event**, for training materials **TrainingMaterial**. Fig. 3 displays a schematic overview of the protocol steps.

![overview of the steps for automated harvesting](img/cpz1682-fig-0003-m.jpg)

*********************

            {{1}}
*********************

1.	Select the profile [**Learning resource**](https://bioschemas.org/profiles/TrainingMaterial/1.0-RELEASE) for the automated harvesting of a specific training material. 
For training events, select **Course/CourseInstance** for training by using the profile described at https://bioschemas.org/profiles/CourseInstance.

In our example, the material is provided by the [NanoCommons initiative](https://nanocommons.github.io/user-handbook/) located [here](https://nanocommons.github.io/tutorials/enteringData/index.html).

2.	Manually annotate the HTML web site containing the training material with JSON-LD markup. 

In our example, we use **index.html**.In order to inspect this object, browse to the [index page](https://nanocommons.github.io/tutorials/enteringData/index.html) of the NanoCommons initiative.  More specifically, enter a JSON object as shown below in a script element with the attribute type="application/ld+json". This script element needs to be placed as a child of the head element of the HTML page. This results in the following script HTML element.

```json
 <script type="application/ld+json">
  {
	"@context": "https://schema.org/",
	"@type": "LearningResource",
	"http://purl.org/dc/terms/conformsTo": { "@type": "CreativeWork", "@id": "https://bioschemas.org/profiles/TrainingMaterial/1.0-RELEASE" },
	"name": "Adding nanomaterial data",
	"version": "0.9.3",
	"description": "This tutorial describes how nanomaterial data can be added to an eNanoMapper server using a RDF format.",
	"license": "https://creativecommons.org/licenses/by/4.0/",
	"keywords": "ontologies, enanomapper, RDF",
	"url": "https://nanocommons.github.io/tutorials/enteringData/",
	"provider": {
  	"@type": "Organization",
  	"name": "NanoCommons",
  	"url": "https://www.nanocommons.eu/"
	},
	"audience": {
  	"@type": "EducationalAudience",
  	"educationalRole": "Graduates"
	},
	"inLanguage": {
  	"@type": "Language",
  	"name": "English",
  	"alternateName": "en"
	},
	"author": [
  	{
    	"@context": "https://schema.org",
    	"@type": "Person",
    	"name": "Egon Willighagen",
    	"identifier": "https://orcid.org/0000-0001-7542-0286",
    	"orcid": "https://orcid.org/0000-0001-7542-0286"
  	}
	]
  }
</script>
```
 
3.	Validate the individual page with the [schema.org validator](https://validator.schema.org) by pasting the URL into the Fetch URL tab. The validation procedure will indicate if you have used non-existing properties of the Bioschemas profile. If error messages are returned, have a look at the [troubleshooting](#12) section below.

Example: https://nanocommons.github.io/tutorials/enteringData/index.html

![screenshot for login page](img/b369eIQ.png)

4.	Create a sitemap listing the material page URL and save it as sitemap.xml. This sitemap.xml file needs to be publicly browsable on the internet. In our example, the manually created sitemap.xml is [published](https://nanocommons.github.io/sitemap.xml).
More complex mechanisms for sitemap creation are available by using content management systems like Drupal and using specific sitemap plugins.

5.	Register a Content Provider in TeSS  following the basic protocol 6. Use ‘NanoCommons’ as a provider.

6.	Make note of your Content Provider's exact Title and URL referring to the properties Title and URL mentioned in basic protocol 6. In our example, it is ‘NanoCommons’ as Title and ‘https://www.nanocommons.eu/’ as URL. 

7.	Go to https://github.com/ElixirTeSS/bioschemas_sources/edit/main/sources.yml to edit the sources.yml file and add your content provider details along with the URL to your sitemap, e.g.

  - title: NanoCommons
      url: https://www.nanocommons.eu/ 
      source: https://nanocommons.github.io/sitemap.xml
 
8.	Commit your change, click the button to open a pull request, which will then be reviewed. After review, the new training material will appear in TeSS. 
Your content should appear in TeSS the following day. Each source will be scraped once per day at ~3AM UTC.

*****************

## Curation of manual content for automatically harvested event or training materials 

In protocols [5](#protocol-5-manual-registration-of-a-training-event-in-tess) or [6](#protocol-6-manual-registration-of-training-material-in-tess) about manual provisioning of the information, it is important to note that much of the content in TeSS is retrieved and kept up-to-date via automated scrapers (see [protocol 8 Registering events/materials automatically](#protocol-8-automated-harvesting-of-a-training-event-in-tess)) that pull information from public web pages and APIs in regular intervals. To prevent TeSS from overwriting a field you have just changed, click the “locker” icon to lock the field. Fields marked as locked will not be overwritten by an automatic procedure.

For the training event and material, several mandatory fields have to be filled in, marked with the “*” symbol. There are also a number of optional fields. 

## Recommendations for automated harvesting of a training event in TeSS

TODO: check on exercises

In step 2 of [Protocol 8](#protocol-8-automated-harvesting-of-a-training-event-in-tess), we recommend to follow the guidelines from the BioSchemas community on strategies on how to markup your internet sites by using this tutorial https://bioschemas.org/tutorials/howto/howto_add_markup. If your website is hosted on github, you can follow https://bioschemas.org/tutorials/howto/howto_add_github instead.

In step 3, in the case that you get an error during the validation procedure by the validator from schema.org, consult the properties definitions in the respective profile on schema.org and bioschemas.org. 
In step 7, you need to have an account on github.com to edit the sources.yml file. It is important to make sure the title and url exactly match the Content Provider's title and URL on TeSS. You will make changes in the github project ElixirTeSS/bioschemas_sources, for which you don’t have write access to. Submitting a change will write it to a new branch in your fork <your github user name>/bioschemas_sources, so that you can send a pull request afterwards.
If you have implemented Bioschemas markup on your website and would like your content to appear in TeSS, see our "Bioschemas sources" repository for detailed information on how to proceed: https://github.com/ElixirTeSS/bioschemas_sources#readme

TODO: add new feature on Provider page!

## Troubleshooting

In case, you experience issues with missing concepts in the EDAM ontology, potentially new properties in the form on the TeSS website, a missing hosting institution or synchronizations issues with the automated procedure, please refer to entries in the table below for troubleshooting.

Sources and Solutions to Potential Errors

| Problem                                                                           | Possible Cause                                               | Solution                                                                                                            |
| :-------------------------------------------------------------------------------- | :----------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------ |
| Missing concept in the EDAM ontology.                                             | The EDAM ontology does not cover all scientific topics.      | Select a term which is closely related or submit a request to add new concepts in the EDAM ontology if appropriate. |
| There are new properties in the TeSS interface for events or material.            | Sometimes new properties are added to the interface in TeSS. | Browse to the documentation on the TeSS web site.                                                                   |
| Your host institution is missing when you are registering an event or a material. | The dropdown of the host institutions is a fixed list.       | Enter a new host organisation in the field.                                                                         |                                 
| By using the automated harvesting protocol, your training material does not appear in TeSS on the next day. | The title and/or URL did not exactly match the Content Provider's title and URL on TeSS  | Double check that the title and url exactly match the Content Provider's title and URL on TeSS.| 
| After using the automated harvesting protocol, your training material does not appear in TeSS the following day | You might have provided a validated Bioschemas markup but the harvesting procedure in TeSS did not recognise all the properties correctly.| Submit an issue on the ELIXIRTeSS/Bioschemas_sources github repository.                                                                                                     |


## Relevant internet resources

http://edamontology.org/: This is the main information site about the EDAM ontology, the ontology of bioscientific data analysis and data management. It is a comprehensive ontology of well-established, familiar concepts that are prevalent within bioscientific data analysis and data management (including computational biology, bioinformatics, and bioimage informatics). 

https://schema.org/: Schema.org is a collaborative, community activity with a mission to create, maintain, and promote schemas for structured data on the Internet, on web pages, in email messages, and beyond. On this web site, you can find 

https://bioschemas.org/: The homepage of the Bioschemas community. Bioschemas aims to improve the Findability on the Web of life sciences resources such as datasets, software, and training materials. Here you can find all the information about the three Bioschemas Profiles implemented in this protocol at https://bioschemas.org/profiles as well as various tutorials at https://bioschemas.org/tutorials.

https://validator.schema.org: This is the site where you can easily validate web sites which have been marked up with structured data as e.g. a JSON-LD object.

https://github.com/ElixirTeSS/bioschemas_sources#readme: This is the README file of the github repository bioschemas_sources for the TeSS portal hosted in the ELIXIR TeSS github organisation. 

https://github.com/TrainTheTrainer/ELIXIR-EXCELERATE-TtT:  This is the github repository ELIXIR-EXCELERATE-TtT hosted in the TrainTheTrainer github organisation. 

https://elixirtess.github.io/TeSS_widgets. This is the github repository for TeSS widgets where several implementation scenarios are documented.


## Publications

Beard, N., Bacall, F., Nenadic, A., Thurston, M., Goble, C. A., Sansone, S.-A., & Attwood, T. K. (2020). TeSS: A platform for discovering life-science training opportunities. Bioinformatics, 36(10), 3290–3291. https://doi.org/10.1093/bioinformatics/btaa047

Garcia, L., Batut, B., Burke, M. L., Kuzak, M., Psomopoulos, F., Arcila, R., Attwood, T. K., Beard, N., Carvalho-Silva, D., Dimopoulos, A. C., Angel, V. D. del, Dumontier, M., Gurwitz, K. T., Krause, R., McQuilton, P., Pera, L. L., Morgan, S. L., Rauste, P., Via, A., … Palagi, P. M. (2020). Ten simple rules for making training materials FAIR. PLOS Computational Biology, 16(5), e1007854. https://doi.org/10.1371/journal.pcbi.1007854

Castro, L., Palagi, P. M., Beard, N., Bioschemas Training Profiles Group Members, ELIXIR FAIR Training Focus Group, The GOBLET Foundation, Attwood, T., Brazas, M. (2022). Bioschemas Training Profiles: A set of specifications for standardizing training information to facilitate the discovery of training programs and resources, preprint https://doi.org/10.1101/2022.11.24.516513

Gray AJ, Goble C, Jiménez RC, The Bioschemas Community Bioschemas: from potato salad to protein annotation. International Semantic Web Conference; Berlin. 2017. Accessed 2019 Mar 11. Available from: https://pdfs.semanticscholar.org/74ec/a9c89622bff731b21b03acb4f2400a0f00fa.pdf

Guha R.V., Dan Brickley, Steve MacBeth (2015) Schema.org: Evolution of Structured Data on the Web: Big data makes common schemas even more necessary. Queue, 13(9), 10–37. https://doi.org/10.1145/2857274.2857276

Ison, J., Kalas, M., Jonassen, I., Bolser, D., Uludag, M., McWilliam, H., Malone, J., Lopez, R., Pettifer, S., & Rice, P. (2013). EDAM: An ontology of bioinformatics operations, types of data and identifiers, topics and formats. Bioinformatics, 29(10), 1325–1332. https://doi.org/10.1093/bioinformatics/btt113

Jupp S. et al. (2015) A new Ontology Lookup Service at EMBL-EBI. In: Malone, J. et al. (eds.) Proceedings of SWAT4LS International Conference 2015
the FAIRsharing Community, Sansone, S.-A., McQuilton, P., Rocca-Serra, P., Gonzalez-Beltran, A., Izzo, M., Lister, A. L., & Thurston, M. (2019). FAIRsharing as a community approach to standards, repositories and policies. Nature Biotechnology, 37(4), 358–367. https://doi.org/10.1038/s41587-019-0080-8

Via, A., Palagi, P.M., Lindvall, J.M., Tractenberg, R.E., Attwood, T.K., The GOBLET Foundation et al. (2020). Course design: Considerations for trainers – a Professional Guide. F1000Res. https://doi.org/10.7490/f1000research.1118395.1

Whetzel, P. L., Noy, N. F., Shah, N. H., Alexander, P. R., Nyulas, C., Tudorache, T., & Musen, M. A. (2011). BioPortal: Enhanced functionality via new Web services from the National Center for Biomedical Ontology to access and use ontologies in software applications. Nucleic Acids Research, 39(suppl), W541–W545. https://doi.org/10.1093/nar/gkr469
