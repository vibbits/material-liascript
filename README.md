<!--
author:   Alexander Botzki
email:    Alexander.Botzki@vib.be
version:  0.1.1
language: en
narrator: US English Female

comment:  DMP Introduction

logo: img/liascript_img/Logo.png

link:     https://cdnjs.cloudflare.com/ajax/libs/animate.css/3.7.2/animate.min.css
link:     https://raw.githubusercontent.com/vibbits/material-liascript/master/img/org.css

debug: true

-->

# Protein Structure Analysis

[Slides](https://liascript.github.io/course/?https://raw.githubusercontent.com/vibbits/material-liascript/master/slides-PSA.md)

# How to write a Data Management Plan (DMP)

Welcome to this self-study course about research data management. In this course you will learn more about how you can manage your research data.

A lot of material of this version of the **Learn to write your DMP** course has been developed by RDM Support of University of Utrecht.

It has been transformed to this open version using [LiaScript](http://LiaScript.github.io) by [VIB Bioinformatics Core](https://www.bits.vib.be), [ELIXIR Belgium](https://www.elixir-belgium.org) and [Helis Academy](https://helisacademy.com/en).

For this online training we have made use of existing (training) material of

* [MANTRA](http://mantra.edina.ac.uk/), RDM training by the University of Edinburgh
* [Managing Data @ Melbourne](http://library.unimelb.edu.au/Digital-Scholarship/training_and_outreach/data)
* [Essentials for Data Support](http://datasupport.researchdata.nl/) by Research Data Netherlands
* [RDM website of University of Amsterdam](http://rdm.uva.nl/)
* [RDM website of Erasmus University of Rotterdam](https://www.eur.nl/researchmatters/research_data_management/)
* [Research Data Services](https://www.tudelft.nl/en/library/current-topics/research-data-management/) of Delft University of Technology
* [Data Management Support Hub](https://www.wur.nl/en/Expertise-Services/Data-Management-Support-Hub.htm) at Wageningen University & Research
* [Digital Curation Centre](http://www.dcc.ac.uk/) (DCC)
* [UK Data Archive](http://data-archive.ac.uk/)
* [Australian Networked Data Services (ANDS)](http://www.ands.org.au/working-with-data/data-management)
* [ORION e-learning course](https://www.orion-openscience.org/news/201912/new-launch-orion-mooc-open-science-life-sciences)

With all questions about this course you can contact: info.rdm@vib.be.

This content is licensed with the CC:BY license.
Please click [here](https://creativecommons.org/licenses/by/4.0/) for more information about the CC:BY license.

---

> **1.** Use the arrow buttons above for navigation
>
> **2.** Turn on your sound for the text output.

__By VIB Bioinformatics Core, ELIXIR Belgium and Helis Academy__

* GitHub: https://github.com/vibbits/material-liascript
* LiaScript: https://liascript.github.io/course/?https://raw.githubusercontent.com/vibbits/material-liascript/master/README.md#1

## Introduction

### Welcome to the course

The course consists of 6 chapters, divided in three categories.
- Prepare
  - Data collection
  - Data documentation
- Handle
  - Data storage
  - Data security
- Share
  - Data selection and preservation
  - Data availability for reuse

Each chapter starts with an introduction and ends with an assignment to write that part of your data management plan that corresponds with what you have just learned. You are currently in the introduction chapter. In this chapter you will learn more about the course and the learning environment. The course ends with chapter 7, 'Rounding up".


**Data Management Plans**

The assignment throughout the course is to fill your own data management plan. At the end of each chapter you will be asked to log into DMPonline. With the content in this course, you should be able to apply this to your research project.

**Questions about the course**

If you have technical questions, please contact bits@vib.be.
If you have content related questions, please contact RDM Support: info.rdm@vib.be. All feedback is welcome, as this is still a beta version. Based on feedback from users, more content may be added or existing content may be changed to a different form.

**Technical requirements**

Some activities use HTML5. Make sure your browser has installed the latest updates. If an activity doesn't work, we recommend you use another browser.

**Licenses and credits**

We wish you a lot of fun with the course and we hope it turns out to be a useful learning experience.

The content of the course is adapted from an online course of [University of Utrecht](https://lll-platform-uu.nl/).

Please note: the content in this course is available under the CC:BY license. Some contents have been adopted from other courses (see [first slide](#1)).



### Why manage your research data?

                --{{0}}--
In this video Katarzyna Biernacka explains what data in a research context is.

                  {{0-1}}
*************************

!?[Why manage data?](https://www.youtube.com/watch?v=XCckz_4mlhU)

CC-BY-4.0: Katarzyna Biernacka, HU Berlin & [Discipline Workshops 2019](http://www.discipline-workshops.com/)

************************

                --{{1}}--
Managing your data effectively is crucial to the success of your research. This doesn't only apply to the immediate context of your thesis or publications. Managing your data is a practice that will benefit you throughout your research career. The following list gives an overview of what benefits are evident.

                  {{1}}
1. **Access, Re-use & Recognition**
   * Facilitating future research by allowing others to build on or add to your research data.
   * Increased citations of research data and of publications based on that data.
2. **Efficiency**
   * Increasing your research efficiency by saving time and resources.
   * Preventing duplication of effort by enabling others to use your data.
3. **Quality & Security**
   * Ensuring the integrity and reproducibility of your research.
   * Ensuring that research data and records are accurate, complete, authentic and reliable.
   * Enhancing data security and minimising the risk of data loss.
4. **Compliance**
   * Meeting legal obligations, restrictions and codes of conduct.
   * Meeting the University policy for research data requirements.
   * Meeting funding body grant requirements.
   * Meeting publisher requirements for data access.


### A case to consider

        {{0-1}}
Marleen is an early career researcher. She completed her PhD about four years ago and is now a postdoctoral research fellow at a different university. Since she obtained her PhD, she has published a number of journal articles based on her doctoral research. Her papers have been cited widely in the literature of her field. But just recently a fellow researcher has questioned her findings. He has gone so far as to suggest that the data on which her research was based is inaccurate. One implication is that the data could even have been falsified. Marleen is confident that her research is valid and that her data is accurate.


    --{{1}}--

    - What steps could Marleen take to verify her research findings?
    - What evidence would she need to demonstrate that she hasn't falsified her data?

       {{1}}
Think about your own research. If someone accused you of research misconduct, would you be in a position to defend your research and reputation? List some strategies you could implement right now that would assist you, should you ever find yourself in Marleen’s situation.

### Data disasters – postcards from the edge

The following are real examples where researchers or data centers have lost crucial data. Could any of these ever happen to you? With good planning you could avoid or reduce the impact of such occurrences.

<iframe src="https://h5p.org/h5p/embed/581076" width="699" height="571" frameborder="0" allowfullscreen="allowfullscreen"></iframe>
<script src="https://h5p.org/sites/all/modules/h5p/library/js/h5p-resizer.js" charset="UTF-8"></script>

### University policy framework for research data

For the Flemish universities, it is important that all researchers honour scientific standards, including the meticulous and ethical treatment of research data.
This policy is intended to set out parameters to safeguard the quality, availability and accessibility of research data within any Flemish university. It provides a basis for evaluating compliance with laws, regulations and codes of conduct. The policy also clarifies the various roles and responsibilities of university staff in managing research data.

The highlights of the policy are:
* Archive (relevant and valuable) research data for a minimum of ten years;
* Store data in a structure that is suitable for long-term preservation and later consultation;
* Provide metadata to describe the data with sufficient clarity to ensure they are findable for further research;
* Make archived research data available for access and reuse at and outside VIB insofar as is reasonably possible;
* Each individual researcher / research leader is responsible to draw up a Data Management Plan (DMP) at the start of the research project and to follow up the agreements made in this plan;
* Scientific directors are responsible for the implementation and monitoring of the University policy framework and for drawing up additional faculty guidelines to this end if needed.

TODO: add link to policy

### Policy in Practise

In this short video Prof. dr. Chantal Kemner explains the importance of good data management for Utrecht University. Chantal is full professor of Biological Developmental Psychology in Utrecht at the faculty of social sciences and since 2013 also at the UMCU.

!?[Data Management at UU](https://youtu.be/f48l4Uca9nA)

### Funder requirements

More and more research funders explicitly require you to consider the management and publication of your research data, both during and after your research project. The European Commission and the Flemish funders FWO have explicit policies on research data management.

**European Commission - Horizon 2020**

The European Commission wants “Horizon 2020 beneficiaries to make their research data findable, accessible, interoperable and reusable (FAIR), to ensure it is soundly managed. Good research data management is not a goal in itself, but rather the key conduit leading to knowledge discovery and innovation, and to subsequent data and knowledge integration and reuse.” Horizon 2020 is the biggest research and innovation program of the European Commission.

[![European Commission - Horizon 2020](img/O_Funders_Screenshot_H2020guidelines.JPG "European Commission - Horizon 2020")](https://www.nwo.nl/en/policies/open+science/data+management)

**FWO**

FWO states that “FWO has made data management a key element of its policy for all support channels provided by the FWO. The FWO expects researchers to pay due attention to this dimension before, during and for at least five years after their research.”

[FWO Overview Data Management Plan](https://www.fwo.be/en/the-fwo/organisation/data-management-plan/)

### Funder guidelines and templates

Most funders require you to write a Data Management Plan. A DMP outlines all key aspects of collecting, storing and managing research data during and after a project. For this they provide you with guidelines, forms, templates and examples. For more information you can download the documents under Resources or check out the websites. You can also contact your faculty Research Support Office:

- [EC – Horizon 2020: guidelines](https://ec.europa.eu/research/openscience/index.cfm)
- [FWO template]()

### Writing a data management plan

By now it should be clear that data needs to be properly managed throughout its lifecycle. The most effective way to do this is to create a Data Management Plan (DMP). This will take into account all the stages of the research data lifecycle. As outlined earlier, each individual researcher or research leader is responsible to draw up a data management plan. He or she should do this at the start of the research project. And during the research you should actively follow up on the agreements made in this plan.

Think about our early career researcher Sasha (introduced in ‘Why manage your research materials and data?’) who needs to defend herself against accusations of researcher misconduct. As well as defending against misconduct accusations, some additional benefits of creating a data management plan include:

- Accessing your data more easily;
- Prioritising and balancing activities relating to research data collection and storage;
- Mitigating data loss;
- Reaching agreement between stakeholders about ownership of data;
- Reducing time and effort in the long term.
The good news is that this online training will take you through the necessary steps to create a plan during the subsequent modules.

### Getting started with DMPonline

We offer you DMPonline to create your Data Management Plan. DMPonline is an international online service that guides you in creating a DMP by answering a series of questions about your research project. It allows you to create, share, store, and revise your data management plans online. You will be asked to complete different sections of your DMP as we go through the other modules. As a result you will have written your own data management plan at the end of this course.

With DMPonline you can:

* Write your plan and keep it up-to-date
  * You can easily update your DMP throughout the lifecycle of a project

* Share plans online
  * DMPonline allows collaborative access, so you can share your DMP with other researchers, within and outside of Utrecht University.

* Create multiple plans
  * You can store different DMPs for different projects. And you can make a copy of a previous plan as the basis for writing a new one.

* Download plans
  * You can download your DMP in a variety of formats.

We recommend that graduate researchers share their data management plans with their supervisor(s).

!?[DMPonline introduction](https://player.vimeo.com/video/251506151)

### About RDM Support

RDM Support provides all kinds of research data management assistance to researchers of VIB in all stages of their research. This can range from one-off individual advice to large-scale infrastructure coordination.

TODO: add links Flemish Universities

## Prepare: Data Collection

### Introduction to data collection

By now you will have obtained some idea of what research data management is all about. Now we will have a more in-depth look into the different phases of your research by starting with data collection.

Data collection involves understanding the different types of data you collect. Depending on the nature of your research, there are different methods of collecting data and thus different types of data.

Your data may be physical (paper records or archival forms) or digital (database contents or Excel data). The source of your data may be external, you collect it yourself or you generate it from a machine.

When you write your data management plan you will need to take into account the type of data you collect, the source of the data, and how you will process and analyse your data.

**In this part of the course you will learn to:**

* identify preferred file formats for your research data;
transform your files into a preferred format;
discover sources for existing data;

* discover existing data yourself;

* how to assess the usefulness of existing data;
understand how a workflow leads to data products;

* estimate costs involved with managing your data;

* check the current and expected costs for your research data;

* write the data collection section for your data management plan.

You can watch the video below, provided by TU Delft, about data collection. The video stops at 1:12.

!?[Data collection](https://www.youtube.com/watch?v=AqnVrnVdv2Y)

### Preferred formats for your research data

                  {{0-1}}
****************************** ![Introduction ](img/00_Preferred-formats.png)
******************************

                 --{{1}}--
This e-module is based on the online Research Data Management training 'MANTRA' of The University of Edinburgh and Managing Data @ Melbourne

                 {{1-2}}
******************************

![Introduction ](img/01_Preferred-formats_Learning_Objective.png)

CC BY: [https://mantra.edina.ac.uk/](https://mantra.edina.ac.uk/)

*****************************

              --{{2}}--
The file formats you use to generate your research data will influence how you can manage them over time, i.e. a program or application must be able to recognise the file format in order to access your data within the file.
For example, a web browser is able to process and display a file in the HTML file format so that it appears as a web page. If the browser encounters another file type, it may need to call on a special plug-in to view it. Or it may simply let you download the file to view if it can recognise it in another program.

                 {{2-3}}
******************************

To identify the file format, files usually have a file name extension, or suffix that follows a full stop in the file name and contains three or four letters, like for example:

TODO: add PDF with links to preferred file formats
* .txt    text
* .pdf    portable document format
* .jpg    joint photographic experts group
* .csv    comma separated values
* .html   hypertext markup language
* .xml  extensible markup language  
* .rtf  rich text format

******************************

               {{3-4}}
******************************

![Proprietary formats](img/02_Preferred-formats_proprietary-formats-01.png)

![Proprietary formats](img/02_Preferred-formats_proprietary-formats-02.png)

******************************

              {{4-5}}    
******************************

Determine which format is proprietary and which is an open format

[[.xml] [.pdf] [.psd] [.odf] [.ppt] [.docx] [.csv] [.xls]]
[[ ] [ ] [X] [ ] [X] [X] [ ] [X]]  proprietary format
[(X) (X) ( ) (X) ( ) ( ) (X) ( )]  open format

******************************

              {{5-6}}    
******************************

TODO: list of preferred formats

******************************

              {{6-7}}    
******************************

![Risks of file conversion](img/03_Preferred-formats-file-conversion.png)

******************************

             {{7-8}}    
******************************

While file conversion or migration sometimes has to be done, there are also risks.

Which ones can you think of?

[[X]] file size may change and even become surprisingly large
[[X]] blanks used as missing data code
[[X]] special characters and end of line returns may change
[[X]] relation among items in a table and among tables may be lost
[[X]] layers, color fidelity and resolution may be lost or changed in image files
[[X]] fonts, footnotes and links to other documents may change
[[X]] frame rate, sound quality, codecs and wrappers may be altered in multimedia files
[[X]] last characters in rows (due to row size limitations) may be altered

******************************

              {{8-9}}    
******************************

Open the following .docx file to the preferred format .txt: [PreferredFormatsExcersizePenguinDOC.docx](img/PreferredFormatsExcersizePenguinDOC.docx)

1. Convert this docx file to the preferred format .text
2. Open the text file in an editor
3. Is all formatting perserved OK?

[( )] Yes
[(X)] No
***********************************************************************

![image](https://upload.wikimedia.org/wikipedia/commons/d/d0/Creative-Tail-Animal-lion.svg)

***********************************************************************

******************************

          {{9-10}}    
******************************

Open the following .docx file to the preferred format .txt: [PreferredFormatsExcersizePenguinDOC.docx](img/PreferredFormatsExcersizePenguinDOC.docx)

1. Convert this docx file to the preferred format .odt
2. Open the .odt file
3. Is all formatting perserved OK?

[( )] Yes
[(X)] No

******************************

         {{10}}    
******************************

![Data compression](img/04_Preferred-formats-data-compression.png)

******************************

### Discovering existing data

TODO: block

### Describe what kind of data you will generate

Having a clear view of what data you will generate will enable you to plan its management. You can create an overview of the data you produce or collect by drawing the data in a workflow, or noting down in a table.

Please watch the video below. Tessa Pronk will explain to you how to describe your data.

!?[Describing your data](https://www.youtube.com/watch?v=KE2UpZY4wYA&feature=youtu.be)

### Order elements in your data flow

TODO: add H5P quiz

### Copyright and Intellectual Property Rights (IPR) issues

Copyright is a form of intellectual property right which arises automatically if an original work is created. Copyright may affect the way data may be stored, shared and reused. You should ask yourself who the copyright holder of your datasets is, especially when you use existing data or when you collaborate with external parties.

**Using someone else’s research data**
SURF provides a brief guide to determining what consent is necessary to reuse someone else’s data (see "A brief guide ... someone else's data" in the resources below)  

**Clarifying the ownership of your research data**
TODO: change accordingly for VIB
Officially VIB, as your employer, is considered the rights holder to the research data you create. You, as a researcher, have the primary responsibility for taking care of the data. Questions on data exploitation may be even more important than those of ownership. Who can use the data? Who can publish it? Who can provide it to third parties?  

We strongly recommend that you deal with the issues around data exploitation at an early stage of your research project. Write down agreements between yourself, your supervisor, project members and other interested parties in your Data Management Plan.

TODO: change accordingly
RDM Support offers you a Guide to legal instruments and agreements for research data management (see the Guide 'Legal instruments and agreements')

**Confidential or privacy-sensitive data**
When your research project has received data under confidentiality or under legal privacy restrictions, you will have to identify and explain how you will deal with these restrictions in your data management plan (also see ‘Learning Unit: Handle - Data security’).

### Costs involved with managing your data

Check
https://www.uu.nl/en/research/research-data-management/guides/costs-of-data-management

        --{{1}}--
The costs of data management and sharing activities must be included into your research, in terms of time and resources needed.

         {{1}}
*****************************

**1. Data Management Cost Guide**

When you plan your research you may not be able to oversee all costs involved. Nevertheless, it is useful to have an idea of possible costs at an early stage. You can use the Guide 'Costs of Data Management', which is a practical overview of possible costs per activity within each phase of the research process. Note: The Cost Guide offers cost indications and examples. These are not real prices.

**2. Budget your data management costs**

You are advised to budget the data management costs as separate data management costs. These costs are eligible for funding with funders like NWO and the European Commission, as long as the costs are invoiced before the end of the project.

**3. Planning can save time and money**

Planning an early start for certain activities within your research project can lower the costs for data management in the run of your project. You can save time by:

Properly describing your data while collecting it, instead of doing it afterwards
Choosing the right file format so that file conversion afterwards is not necessary
Hiring an experienced data manager
Spending time to think about data activities beforehand can help prevent unexpected extra efforts and costs later on in your research project.

*****************************

### Check the current and expected costs for your research data

You have just learned that in many parts of a research project there are data related costs. These costs depend on the type and volume of data you produce, analyse and store.

check link to file (calculation)
https://lll-platform.uu.nl/pluginfile.php/4907/format_elevated/resource/0/Cost%20overview.docx

### Write your data management plan for your data collection

Go to DMPonline and open your draft data management plan created in the Introduction.

You have now completed the module Data collection. You should be able to complete the following questions in the section Data collection:

* Will you use existing data?
* What data will you collect or create?
* How will the data be collected or created?
* How will you manage rights issues?
* What are the costs involved in managing and storing your data?

## Prepare: Data Documentation

### Introduction to documentation and metadata

By now you understand how to describe your data collection in terms of, for example, type, size, and format. You have identified this for your own research data.

Now we will look into the documentation and metadata which will accompany your data. Documentation and metadata are essential to understand what a dataset means and to make it reusable in the future.

In this part of the course you will learn to:

* Understand why your research data needs to be documented and why other researchers rely on this documentation;
* Understand why metadata is used and when you need to use this;
* List different descriptions of data (record, study, code);
* Describe the purpose/goals of descriptions of data;
* State why standards are important and how to find and apply them;
* Name and apply the do’s and don’ts in arranging folders and naming files;  
* Recognise the importance of managing data files;
* Apply the gained knowledge about data documentation to your own research data and write the data documentation section for your data management plan.

### Introduction to documentation and metadata

          {{0-1}}
****************************** ![Introduction ](img/00_Metadata.png)
******************************

           {{1-2}}
******************************

![Introduction ](img/01_Metadata_Learning_Objective.png)

*****************************

            --{{2}}--
Tips for data documentation - John MacInnes, professor of Sociology of the University of Edinburgh, explains why it is necessary to document each step of your research and how this will benefit you in the long term.

              {{2-3}}
******************************

!?[John MacInnes Tips on Documentation](https://youtu.be/EIZsxT-fIiQ)

******************************

               {{3-4}}
******************************

**Examples of data documentation**

Since there is a wide variety of types of data and types of research, there are many different ways of documenting data. A few examples of data documentation are:

* Laboratory notebooks and experimental procedures
* Questionnaires, codebooks, data dictionaries
* Software syntax and outout files;
* Information about equipment settings & instrument calibrations
* Database schemes
* Methodology reports
* Provenance information about sources of derived or digitised data

******************************

              {{4-5}}
******************************

What data documentation will you use?
TODO: quiz without an answer

[[]] Laboratory notebooks
[[]] Code on github
[[]] Software syntax or outout files
[[]] Questionnaires
[[]] Experimental protocols
[[]] Information about equipment settings & instrument calibrations
[[]] Database schemes
[[]] Methodology reports
[[]] Provenance information about sources of derived or digitised data

******************************

          --{{5}}--
There are many different ways to set up and organise your documentation.

          {{5-8}}
*****************************

**Project level**

Project level documentation documents what the study sets out to do; how it contributes to new knowledge in the field, what research questions/hypotheses are, what methodologies are used, what samples are used, what intruments and measures are used, etc. A complete academic thesis normally contains this information in details, but a published article may not. If a dataset is shared, a detailed technical report needs to be included for the user to understand how the data were collected and processed. You should also provide a sample bibliographic citation to indicate how you would like secondary users of your data to cite it in any publication.

******************************

          {{6-8}}    
******************************

**File or database level**

File or database level documentation documents how all the files (or tables in a database) that make up the dataset relate to each other, what format they are in, whether they supersede or are superseded by previous files, etc. A readme.txt file is the classic way of accounting for all the files and folders in a project.

******************************

           {{7-8}}    
******************************

**Variable or item level**

Variable or item level documentation documents how an object of analysis came about. For example, it does not just document a variable name at the top of a spreadsheet file, but also the full label explaining the meaning of that variable in terms of how it was operationalised.

******************************

            --{{7}}--
John MacInnes, professor of Sociology of the University of Edinburgh, speaks about how data documentation can help to find a way in often voluminous data collections of different copies, routings, syntaxes, samplings, etc.

           {{8-9}}    
******************************

**On the necessity of data documentation in secondary data analysis**

!?[John MacInnes, Data documentation in secondary analysis](https://youtu.be/Ebaiwg08CW8)

******************************

            {{9-10}}
************************

Looking back at your previous research project: Did you ever have problems reusing other people's data because of lack of documentation?

[[]] Never tried
[[]] Successfully reused
[[]] Had to ask clarification
[[]] Had to abandon the reuse attempt

****************************

             {{10-11}}
****************************

![Lab notebooks](img/02_Metadata_Lab-Notebook.png)

Thorough and effective management of laboratory data and the routine documentation of all lab procedures is a highly important responsibility for all researchers.

TODO: Link to eLab Notebook

*****************************

            {{11-12}}
************************

**An introduction to metadata**

Watch this web lecture to learn about the different types of metadata and how metadata can help make your research data better findable. You are pointed to useful sources for metadata standards.

!?[Ins and outs of metadata and data documentation](https://youtu.be/h0oZ3swbTJ0)

****************************

          {{12-13}}
************************

**identify different types of metadata**

HP5 quiz or matrix quiz

****************************

           {{13-14}}
************************

**Metadata for different disciplines**

Different disciplines like biology, earth sciences, physical sciences and social sciences and humanities have their own standards. By choosing a well-supported standard, you will maximise the chance that your data can be re)used and understood by other researchers.

!?[Metadata standards](https://youtu.be/AvL7hEk8RJQ)

****************************

       {{14-15}}
************************

**Metadata for different disciplines**

Useful links to metadata standards:

* [Biology](http://www.dcc.ac.uk/resources/subject-areas/biology)
* [General Sciences](http://www.dcc.ac.uk/resources/subject-areas/general-research-data)

A community-maintained [directory of metadata schemas](http://rd-alliance.github.io/metadata-directory/) which has been set up under the auspices of the Research Data Alliance.

A list of metadata standards and other standards developed by [FairSharing](https://fairsharing.org/).

*****************************

           {{15-16}}
************************

**Controlled vocabulary**

![Controlled vocabulary](img/03_Metadata-controlled-vocabulary.png)

****************************

            {{16-19}}    
******************************

**Improve a record description**

Take a look at the record descriptions n the table below and answer the question below and in the following pages.

| Soil Sample       | Condition     | Length| Classx |
| ----------------- |:-------------:| -----:|:-------|
| A1                | low           | $458  | III    |
| A2                | low           | $391  | II     |
| A3                | medium        | $422  | IV     |

x according to the classification from last experiment

Sample - is the value clear?

[(X)] Yes
[( )] No
******************************

                    {{16}}
TODO: explain why

******************************

******************************

                   {{17-19}}    
******************************

**Improve a record description**

Take a look at the record descriptions n the table below and answer the question below and in the following pages.

| Soil Sample       | Condition     | Length| Classx |
| ----------------- |:-------------:| -----:|:-------|
| A1                | low           | $458  | III    |
| A2                | low           | $391  | II     |
| A3                | medium        | $422  | IV     |

x according to the classification from last experiment

Condition - is the value clear?

[( )] Yes
[(X)] No
******************************

        {{17}}
Correct! It is not clear what low or medium as condition means.

******************************

******************************


                   {{18-19}}    
******************************

**Improve a record description**

Take a look at the record descriptions n the table below and answer the question below and in the following pages.

| Soil Sample       | Condition     | Length| Classx |
| ----------------- |:-------------:| -----:|:-------|
| A1                | low           | $458  | III    |
| A2                | low           | $391  | II     |
| A3                | medium        | $422  | IV     |

x according to the classification from last experiment

Length - is the value clear?

[( )] Yes
[(X)] No
******************************

        {{18}}
Correct! It is not clear what is meant by length. Also a unit for the values is missing. Is it meters, centimeters, or seconds?

******************************

******************************

                   {{18-19}}    
******************************

**Improve a record description**

Take a look at the record descriptions n the table below and answer the question below and in the following pages.

| Soil Sample       | Condition     | Length| Classx |
| ----------------- |:-------------:| -----:|:-------|
| A1                | low           | $458  | III    |
| A2                | low           | $391  | II     |
| A3                | medium        | $422  | IV     |

x according to the classification from last experiment

Class - is the value clear?

[( )] Yes
[(X)] No
******************************

        {{18}}
Correct! There is a reference that the classes are explained somewhere. But no link to the document is given.

******************************

******************************

### Data standards explained

Your dataset can be standardised in various aspects. Standardisation, in general, makes data comparable and interpretable. In other words, your data becomes interoperable by applying standards. Datasets can be combined, compared or are simply easier to reuse. You have to plan standardisation, as it is for many aspects hard or impossible to apply afterwards.

Standardise as much as possible between you and your collaborators or research group. If there are standards established and used in your field of research you are advised to use these.

Here is a list of things you can standardise in your research.

* Standardise how, what and when you measure things by standardising your protocol, or methods and materials For instance, is there a standard set of questions for ‘quality of life’? Is there a standard procedure to house mice for your purpose? What aspects do you measure? At what parameter values (age, concentration, etc.)? When do you measure (every two hours, every gram of weight gain, etc.)?

* Standardise your file formats so you can easily exchange results without technical difficulties. Check for standard taxonomies or coding systems within your research discipline.

* Standardise the units in which you note down your results. For instance, do you use mm, cm, m? It is extra work to transform units between experiments.

* Standardise the metadata you use to describe your records or study. What fields will fill in by default, and according to what standard do you define the fields’ names? Will you design a metadata spreadsheet where you specify all things that you will note down?

* Standardise the vocabulary you use. If everyone has the same terminology, it can avoid confusion or misinterpretation. Check for standard taxonomies or coding systems within your research discipline.

### Check your knowledge on standards

Follow the links below for examples of standards. What type of standardisation do the links refer to?

* [Demographic market research](http://www.amplituderesearch.com/market-research-questions.shtml)
* Find via Google: “general morphology score (GMS)”
* [Marine Geoscience Data](http://www.marine-geo.org/submit/guidelines.php)
* [International Union of crystallography](http://www.iucr.org/resources/cif/spec/ancillary/abbreviations)
* [The Cultural Objects Name Authority](http://www.getty.edu/research/tools/vocabularies/cona/index.html))
* [SI Units](https://www.nist.gov/pml/weights-and-measures/metric-si/si-units)
* [UK data service](https://www.ukdataservice.ac.uk/manage-data/format/recommended-formats)

TODO: add H5P exercise

### Folder structure and file naming

          {{0-1}}
******************************

![Introduction ](img/00_Folder-structure.png)

******************************

           {{1-2}}
******************************

![Introduction ](img/01_Folder-structure-Learning-Objective.png)

CC BY: [https://mantra.edina.ac.uk/](https://mantra.edina.ac.uk/)

*****************************

              {{2-3}}
******************************

![Introduction to good file management](img/02_Folder-structrue-introduction-file-management.png)

******************************

              --{{3}}--
Trying to find a data file taht you need which has been stored or named incorrectly or inaccurately can be both frustrating and a waste of valuable time. In this short video Jeff Haywood, professor at the University of Edinburg, explains his experiences with good and bad file management.

               {{3-4}}
******************************

!?[Jeff Haywood Importance of good file management](https://youtu.be/i2jcOJOFUZg)

******************************

          {{4-8}}
*****************************

**Project level**

Project level documentation documents what the study sets out to do; how it contributes to new knowledge in the field, what research questions/hypotheses are, what methodologies are used, what samples are used, what intruments and measures are used, etc. A complete academic thesis normally contains this information in details, but a published article may not. If a dataset is shared, a detailed technical report needs to be included for the user to understand how the data were collected and processed. You should also provide a sample bibliographic citation to indicate how you would like secondary users of your data to cite it in any publication.

******************************

          {{6-8}}    
******************************

**File or database level**

File or database level documentation documents how all the files (or tables in a database) that make up the dataset relate to each other, what format they are in, whether they supersede or are superseded by previous files, etc. A readme.txt file is the classic way of accounting for all the files and folders in a project.

******************************

           {{7-8}}    
******************************

**Variable or item level**

Variable or item level documentation documents how an object of analysis came about. For example, it does not just document a variable name at the top of a spreadsheet file, but also the full label explaining the meaning of that variable in terms of how it was operationalised.

******************************

            {{8-9}}
******************************

**Choose the best chronological file name**

Which of the file names below is the most appropriate?

[[X]] 2019-03-24_Attachment
[[ ]] 24 March 2006 Attachment
[[ ]] 240306attach
******************************
               {{8}}
Correct! Using a date in the format Year-Month-Day will maintain the chronological order of your files.

******************************

******************************

            {{9-10}}
******************************

**Choose the best descriptive file name**

Which of the file names below is the most appropriate file name?

[[ ]] labtox_recent_110810_old_version.sps
[[X]] 2010-08-11_bioasssay_tox_V1.sps
[[ ]] FFTX_3776438656.sps
******************************
               {{9}}
Correct! Keep the file names short and relevant while using sufficient characters to capture information. Do not name files recent or final or definitive_final, a date or version number will suffice.

******************************

******************************

            {{10-11}}
************************

![Batch renaming](img/03_Folder-structure-batch-renaming.png)

****************************

             {{11-12}}
****************************

![Lab notebooks](img/04_Folder-structure-version-control.png)

*****************************

           {{12-13}}
******************************

**How would you treat your data**

Why should you discard or delete obsolete versions of data?

[[ ]] The most current version is the only relevant version.
[[X]] You have several versions of files in a state between versions
[[ ]] You are exceeding the storage space available to you.
******************************
               {{12}}
Correct! Too many similar or related files may be confusing to yourself and to anyone else wanting to access or use your data. You may think that you know which data file is which but that may not always be the case as time passes and the number of different versions increases. It is easier to maintain a manageable number of versions with a clear naming structure. As long as the original raw or definitive copy is retained and processing is well documented, the intermediate working files can and should be discarded.

******************************

******************************

       {{14-15}}
************************

**Fill the blanks**

TODO: add H5P

*****************************

### Write your data management plan for your data documentation

Go to DMPonline and open your draft data management plan created in the Introduction.

You have now completed the module Data documentation. You should be able to complete the following questions in the section Data documentation:

* How will you structure your data?
* How will the data be described and documented?
* What standards will you use?

## Handle: Data storage

### TODO: specific chapter on storage

### Write your data management plan for your data storage

Go to DMPonline and open your draft data management plan created in the Introduction.

You have now completed the module on data storage. You should be able to complete the following questions in the section ‘Data documentation’:

Where will you store your data?
How will the data be backed up?
After finishing this part in DMPonline, please return to the learning environment and click on [Complete]. This takes you back to the course overview. Continue with the next learning unit.

You can ask your faculty or project data manager or RDM Support for a review of your DMP once you have finished writing all or parts of your DMP.

## Handle: Data security

### Introduction to data security

By now you know more about how to manage your data collection, how to organise and document your research data and where and how to store your data.

Now we will take you into the world of keeping data safe and secure.

In this part of the course you will learn to:

* Understand different ways data breaches can happen so you know how to avoid this;
* Understand different ways to achieve data security and reasons for access restrictions;
* Link situations to different legal contracts that arrange access and security;
* Recognise indirectly identifiable data;
* Store and manage privacy-sensitive data;
* Apply the gained knowledge about data security to your own research data and write the data security section for your data management plan.

**Loss of data, loss of academic career**

The loss of scientific data can have a devastating impact on careers. Imagine that you loose all of the research data you've been diligently collecting for four years. Now imagine the knock-on effect: you won't get the PhD you've been working towards, affecting your future career. This nightmare happened to Billy Hinchen, a biologist at Cambridge University. Listen to his story.

!?[Billy Hinchen data loss](https://youtu.be/3xlax_Iin0Y)

### Data breaches

There are several examples of (mainly online) data storage going wrong, leading to leaks of sensitive and personal information.

The picture below shows the biggest cases of data breaches in the past 10 years. They involve some well-known, highly regarded and trusted companies as well as some practices from the academic world.

[![examples about data breaches](img/data-breaches.png)](http://www.informationisbeautiful.net/visualizations/worlds-biggest-data-breaches-hacks/)

### Prevent unauthorised access

Data security may be needed to protect intellectual property rights, commercial interests, or to keep personal or sensitive information safe. Data security involves security of data files, computer system security and physical data security. All three need to be considered to ensure the security of your data files and to prevent unauthorised access, changes, disclosure or even destruction. Data security arrangements need to be proportionate to the nature of the data and the risks involved. Attention to security is also needed when data are to be destroyed.  If data destruction is in order, you need to make sure that the destruction process is irreversible.

Learn about different measures depending on the kind of security you need.

            {{1}}
**Security of data files**

            {{2-3}}
**************************

The information in data files can be protected by:

* Controlling access to restricted materials with encryption. By coding your data, your files will become unreadable to anyone who does not have the correct encryption key. You may code an individual file, but also (part of) a hard disk or USB stick
* Procedural arrangements like imposing non-disclosure agreements for managers or users of confidential data
* Not sending personal or confidential data via email or through File Transfer Protocol (FTP), but rather by transmitting it as encrypted data (eg. [FileSender](https://filesender.belnet.be) )
* Destroying data in a consistent and reliable manner when needed
* Authorisation and authentication: for personal data you have to give very selective access rights to specified individuals.

**************************

            {{3}}
**************************

**Computer security systems**

***************************

            {{4-5}}
***************************

The computer you use to consult, process and store your data, must be secured:

* Use a firewall
* Install anti-virus software
* Install updates for your operating system and software
* Only use secured wireless networks
* Use passwords and do not share them with anyone. Do not use passwords on your UU computer only, but also on your laptop or home computer. If necessary, secure individual files with a password.
* Encrypt your devices (laptop, smartphone, USB stick/disk).

****************************

            {{5}}
**************************

**Physical data security**

***************************

            {{6-7}}
***************************

With a number of simple measures, you can ensure the physical security of your research data:

* Lock your computer when leaving it for just a moment (Windows key + L)
* Lock your door if you are not in your room
* Keep an eye on your laptop
* Transport your USB stick or external hard disk in such a way that you cannot lose it
* Keep non-digital material which should not be seen by others, in a locked cupboard or drawer.

****************************

            {{7}}
**************************

**Data classification**

***************************

            {{8-9}}
***************************

TODO: what to do with classified data

****************************

            {{9}}
**************************

**Data that contain personal information**

***************************

            {{10-11}}
***************************

These data should be treated with higher levels of security than data which do not. You will learn more about privacy-sensitive data in the e-module.

****************************

### What is your experience with unauthorised access to your research data?

We are interested to know if you have ever experienced unauthorized access to any of your research data. When you give your reply, we will show you an overview with the responses of other researchers in this course. All responses will be processed anonymously.

[(1)] No, I am sure about that
[(2)] Not that I am aware of
[(3)] Yes, without much consequences
[(0)] Yes, with severe consequences

### Legal agreements and contracts

Often other people are required to handle your data, or you might be the person that handles other people’s data.

To arrange the security of the research data you work with, in many cases you have to make a (legal) agreement with other people involved. These agreements will make explicit permitted uses, retention time, and agreed upon security measures. Find out what legal contracts you can use by studying the figure below.  TODO: Visit the Guide 'Legal instruments and agreements' for more information

For tailored advice and templates, contact Legal Affairs via your faculty Research Support Officer (RSO)
TODO: add link
![Legal Agreement contacts](img/AgreementsPicture.png)

### When to use which legal contract?

You have been acquainted with the different flavors of legal agreements. Is it clear to you when you need which agreement? Please answer the following questions by choosing the right kind of agreement.

TODO: add quiz or H5P quiz

### Privacy-sensitive data

             {{1-2}}
**************************

![start privacy-sensitive data](img/00_privacy-sensitive-data.png)

**************************

             {{2-3}}
**************************

![start privacy-sensitive data](img/01_privacy-sensitive-data-learning-objectives.png)

***************************

             {{3-4}}
**************************

**Privacy in a nutshell**

---

Privacy is a fundamental right. With regards to privacy, we all have two perspectives:

1. How is your privacy protected?
2. How can we, as a researcher, protect the privacy of the people involved in our research (the data subjects)?

TODO: add link to document and image screenshot
![privacy reference card](img/LCRDM-privacy-reference-card -why-Version-02.pdf)

***************************

           {{4-5}}
***************************

**Six principles from the European General Data Protection Regulation 1/2**

The European General Data Protection Regulation (GDPR) outlines how we should work with privacy-sensitive data.

TODO: create working infographics
see http://gdprcoalition.ie/infographics

**************************

           {{5-6}}
***************************

**Six principles from the European General Data Protection Regulation 2/2**

According to the GDPR  processing of personal data must be done according to 6 principles.

TODO: create HP5 document
see UU learn and Understanding the GDPR University Groningen

**************************

           {{6-7}}
***************************

**Privacy by design**

To comply with the six principles from the GDPR, you can implement privacy by design. This means that you design a data management plan with measures on both IT and procedural level.

TODO: Video privacy by design?

**************************

           {{7-8}}
***************************

**Which data breach is breached?**

Can you recognise the principles that are breached in the different ways personal data is processed?

TODO: H5P quiz 7 cases

**************************

           {{8-9}}
***************************

**Storing personal data 1/2**

![storing personal data](img/02_privacy-sensitive-data-personal-data-01.png)

**************************

           {{9-10}}
***************************

**Storing personal data 2/2**

Only if the access can be unambiguously be restricted to authorised persons, can data be stored without such measures.

Should you want an elaborate visualisation of what is considered identifiable data, check out the information sheet at the Future Privacy Forum.

[visual guide to practical data de-identification](https://fpf.org/2016/04/25/a-visual-guide-to-practical-data-de-identification/)

**************************

           {{10-11}}
***************************

**Can you recognize identifiable data?**

[[X]] a collection of GPS data of daily routines
[[ ]] a list of households sizes associated with number of pets
[[X]] MRI scans without identifying metadata.
[[X]] audio recordings with no metadata and no names of the recorded persons
[[ ]] transcripts of interviews without any directly identifying information
[[ ]] a list of gender and grades for a de-identified course
******************************
               {{10}}
Correct! GPS data holds information on where people go. In a daily routine, the track ends at a particular location which is likely the home of the subject. AN MRI scan from the profile of the head can be identifiable. Audio recordings can be identifiable from the tone of the voice. A list of surnames in itself is not identifying nor personal information.

****************************

****************************

            {{11-12}}
***************************

**Access to privacy-sensitive data**

If and how you can make personal data available, depends n the level of sensitivity of your data. The more sensitive, the more restrictions and safeguards need to be put in place to make sure the data does not fall into the hands of unauthorised persons both during and after research.

To determine where the privacy risks lie for your data you will have to do a Data Privacy Impact Assessment (DPIA).
For more information:
TODO:
link to: https://www.uu.nl/en/research/research-data-management/guides/handling-personal-data

Towards the data subjects, you need to be transparent regarding the possible reuse, or retaining of the data for verification requirements, and get their prior consent.

***************************

            {{12-13}}
***************************

**Cases on how to make personal data accessible**

Case 1: YOUth cohort study

YOUTH COHORT STUDY
YOUth (Youth Of Utrecht) is a large-scale, longitudinal cohort following children in their development from pregnancy until early adulthood.

A total of 6,000 babies and children from Utrecht and its surrounding areas will be included in two different age groups and followed at regular intervals.

The YOUth data enables researchers to look for answers to all sorts of scientific questions on child development. A few examples of YOUth data: human bodily material, hours of videos, MRI images, questionnaires, ultrasounds and IQ scores. YOUth encourages and facilitates data sharing. It is one of the leading human cohorts in FAIR and open data in the Netherlands.

More information at: https://www.uu.nl/en/research/youth-cohort-study

Case 2: TODO: other example from Wings?

***************************

            {{13-14}}
***************************

**An introduction to informed consent**

In the module 'Legal agreements and contracts' you learned about informed consent. Informed consetn is very important when working with data which is in any way related to people.

TODO: add graphics on informed consent

One thing to arrange in your informed consent is the possibility for future use, for verification or reuse. In your informed consent, it is important to be clear on future use of data.

***************************

            {{14-15}}
***************************

**Informed consent for data sharing**

One thing to arrange and to be crystal clear about in your informed consent is the possibility for future use of your data, for verification or reuse.

Check the sentences that do permit data sharing if used as a single statement.

[[X]] Any personal information that reasonably could identify you will be removed or changed before files are shared with other researchers or results are made public.
[[X]] Other genuine researchers (may) have acces to tis data only if they agree to preserve the confidentiality on the information as requested in this form.
[[ ]] Any data that could identify you will be accessible only to the researchers responsible for performing this study.
[[ ]] All personally identifying information collected about you will be destroyed after the study.
******************************
               {{14}}
Correct!

Sharing of research data that relates to people can often be achieved using a combination of obtaining consent, anonymizing data and regulating data access. If the statement towards the data only mentions the current study, sharing is not explicitly possible. You should add some sentence to make it clear to participants that the data could be used for further research, deidentified where possible, or identifiable with enough safeguards and security measures, if it is not.
******************************

******************************



***************************

### Write your data management plan for your data security

Go to DMPonline and open your draft data management plan created in the Introduction.

You have now completed the module on data security. You should be able to complete the following questions in the section ‘Data security’:

* Will you use or collect any confidential or privacy-sensitive data?
* What measures will you take to ensure the security of any confidential or privacy-sensitive data?
* What measures will you take to comply with security requirements and mitigate risks?
To whom will access be granted/restricted?

## Share: Data selection and preservation

### Introduction to data selection and preservation

Research should be transparent and you should always be able to revert back to your data if necessary and be able to show others how you came to your results. Therefore, your research data with all information reasonably necessary for verification needs to be preserved.

With well-managed and preserved research data, you can defend yourself against allegations of mistakes. You can also prevent wrong conclusions from further spreading into the scientific community if there really are mistakes.

In this part of the course you will learn to:

* Select what part of your data should be  preserved for verification purposes;
* Understand the benefits of preserving your data in a public data repository;
* Preserve you data technically correct.

### Long term data preservation

Research data can be preserved for different reasons such as verification and/or possible reuse. It can be your own wish or that of your university, funder or journal.

**Verification**
TODO: adapt this part

The Netherlands Code of Conduct for Academic Practice (VSNU) states that raw data from research must be kept available for a minimum of ten years. This statement is also included in the Utrecht University Policy framework for research data: “Archived research data are to be retained for a minimum of ten years, commencing from the date that the research results are published.”

**Reuse**
It may be worthwhile to make (part of) your data available for a longer period of time and/or for a wider audience. Data which are suitable to keep for reuse are interpretable data on which new research can be based,  independent of the publication.

On the one hand, making research data reusable will need extra effort. On the other hand, possible reuse, even by your future self, might bring you lots of benefits and credits. Consider if your data is worth the effort of making it reusable or if preserving and archiving for verification is enough.

Reuse is explained more in depth in the next part of this course; ‘Availability for reuse’. In this part we will focus on selection and preservation of research data for verification purposes.

### Data package

Keeping data for verification serves the specific goal of having transparent, reproducible research.

**Alternatives to preserving raw data**
If preserving your raw data poses problems, alternatives can also ensure verfication. For instance, transcripts of recorded interviews could hold all important information and may be less privacy-sensitive, so it is reasonable to preserve those instead of the recordings themselves. Also, if raw data is very large, preserving your data only in some processed form could be an alternative. Combined with, for instance, a demonstrable quality check on the processing.

**The contents of your data package**

TODO: add image for illustration/zenodo?

Others should be able to understand what you did. It is not enough to just provide data. Instead you should preserve a package with everything included that is necessary to reproduce your results. Think of including the following:

* Primary (raw) data;
* Secondary (processed) data;
* Protocols;
* Computer code/scripts;
* Lab journals;
* Metadata and/or codebooks describing the data;
* An overview of what the contents of the data package stating what file contains what information, and how these are related.

The data should contain a reference to any publication which is based on the data.

To make understanding your data less dependent on information in the publication, you can also add information on:

* Collection methods;
* Procedures;
* Experimental protocol;
* Your research question;
* Stimuli used;
* Sample descriptions.

This is especially practical if the data package can be found and used on its own account. This is the case if it is published in a data repository or data journal as a data package for reuse.

Do not forget to explicitly state who is responsible for the content of the data package, who is to be contacted in case of a request for access, and under what conditions access is granted.

### Where to preserve what type of data?

During your research, you generate research results that can be made available for others.

A paper or publication is the most traditional way of making results available, but it is by no means the only way. A relatively new way of making results available is using a public data repository.

As you have just learned, preserving your data may serve the purpose of verification or  reuse. Public data repositories cater to both needs. In addition, they handle requests to view or use your data which means you do not have to take care of such requests yourself.

In the example below, you find a workflow for experimental research. What information can be made available in what place? Drag the items on the right to the correct place in the figure. Please note that some items can be used more than once.

TODO: add H5P quiz and PDF solution?

### Accounting for data of others

If you are permitted to use data from other parties, you will have to account for those as well if your research is to be verifiable and reproducible by others. You may recognise this from chapter 1 of this course: Data collection: Discover existing data, weblecture ‘Assessing usefulness of research data of others’ (5 of 10).

You have the following options:

If the used data is preserved correctly somewhere for the coming ten years, refer to the data repository in question;
If it is not taken care of, contact the responsible persons, negotiate correct preservation in a data repository for ten years, and refer to that repository.
If this isn’t possible, try to arrange a local copy that you preserve yourself;
If this isn’t allowed, you will not be able to present the data in case of questions. Therefore, you should question yourself whether you can actually use the data.

![alt-t](img/Cont_5_Share_SelectPreserve_Chart10years.png)

**Accounting for data of others on websites**

If you find interesting information on a website that you want to refer to, it is possible that this information will not be future proof.

The link or web address might change over time (link rot). Or the information on a website is updated, changed or replaced with other content (content drift).

It is possible to archive web pages on a web archive like the [Internet Archive](https://archive.org/web/). You can capture a web page as it appears now for use as a trusted citation in the future (save a page). You will get an alternative link, pointing to the archived, static version of the page. Use this alternative link as a reference to the online information.

### How to preserve your data correctly

In order for the data to survive for the long term, an active preservation regime has to be applied. The bad news is, data automatically gets lost over time.

There are five main ways your data can be lost:

* Digital sources degrade over time ('bit rot');
* File formats and software become outdated;
* The media on which your data is stored becomes outdated or defective;
* Disaster strikes the storage location;
* The person that understands the data finds another job or data simply becomes forgotten.

In this video below you will learn how to minimise the risk of losing data. You are also given good preservation practices.

!?[Preserving data](https://www.youtube.com/watch?v=qENaO0Lk6eo)

### Match the solutions to the data loss

From the weblecture you learned how to prevent data loss. Can you recall all applicable active regimes, as explained in the weblecture?

Below you see a list of solutions to prevent data loss. Underneath that list you see a list of risks for data loss. Please add the number of each solution to the correct risk.

**Solutions to prevent data loss**

1. Have multiple copies. Use a checksum to identify faulty copies
2. Use preferred file formats that can be opened by a wide range of software. Update the file format to a current one.
3. Move data to fresh media well before the media’s expiration date.
4. Have multiple copies. Move data to fresh media well before the media’s expiration date.
5. Document your data well.
6. Advertise the content in a data catalogue.

TODO: add quiz text solution

### Write your data management plan for your data preservation

Go to DMPonline and open your draft data management plan created in the Introduction.

You have now completed the module on data selection and preservation. You should be able to complete the following questions in the section ‘Data selection and preservation’:

* Which data should be preserved and/or shared?
* How and where will you keep your data for the long term?

## Share: Data availability for reuse

### Introduction to data availability for reuse

Thanks to information and communication technology and globalisation new opportunities arise to exchange results of scientific research - publications and research data - and even of scientific methods and practices. This new way of practising science is called ‘open science’.

Open data is a part of this movement towards open science. It is the ambition of universities, governments, funders and publishers to make research data optimally suited for reuse.

There are different reasons why you may not be able to share your research data. Thinking about these issues and challenges when developing your data management plan will help you reflect on such reasons in an early stage.

In this part of the course you will learn to:

* Discern various repositories to store and/or share your data;
* Understand that you should manage a license for use of data yourself;
* Identify ways to ensure that privacy-sensitive data can be shared
* Understand how shared datasets can be cited.

**How frustrating a data request can be**

Not being prepared to share your data can lead to problems in using the data. In this short video, you see what shouldn't happen when a researcher makes a data sharing request! Topics include storage, documentation, and file formats. A made up, yet not unrealistic story.

!?[Data Sharing 3 Short Acts](https://youtu.be/66oNv_DJuPc)

### Introduction to data repositories

In order to preserve, manage, and provide access to your research data, you can deposit your data in a data repository. Data repositories allow permanent access to datasets in a trustworthy environment and enable search, discovery, and reuse of the data they host.

Click on the topics below to find out more about data repositories.

TODO: add repositories from Elixir

                 {{1}}
************************************************

**A wide variety**

************************************************

                 {{2-3}}
************************************************

There is a wide variety of data repositories. Most have the option to publish your dataset using a persistent identifier and some provide the service of long-term preservation. Some repositories host data from various disciplines and others are domain- or discipline specific.

************************************************

                 {{3}}
************************************************

** Choosing a data repository**

************************************************
                 {{4-5}}
************************************************

When choosing a repository for your data be sure to check if the repository meets your criteria or the criteria set by your funder or journal editors.

Criteria to select a certain repository can be:

* Is the repository certified with a [CoreTrustSeal](https://www.coretrustseal.org/) or Data Seal of Approval?
Repositories with a Data Seal of Approval are recognised in the community as a trustworthy source of data.
* Is long term archiving guaranteed or not?
Some repositories will guarantee the legibility of the data, even if the hardware and software become obsolete.
* What are the costs per dataset or gigabyte?
Repositories differ in their cost model, some allow free deposits up to a certain amount of storage
* What is the physical storage location of data?
The location of your data determines under which data protection law it falls. Some repositories store data in the US and others in the EU.
* What is the default license?
Some repositories allow for open or restricted access, or you can specify which license for use you want for your data.

You can use this [repository selection tool](https://www.uu.nl/en/research/research-data-management/tools-services/tools-for-storing-and-managing-data/decision-aid-data-repositories) to help you select a suitable repository.

************************************************

                 {{5}}
************************************************

**Registry of research data repositories**

************************************************

                 {{5-6}}
************************************************

You can browse or search for a data repository in re3data.org. This is a global registry of research data repositories covering different academic disciplines. You can search or browse by subject, content type or country. You can filter the search and browse results on criteria for choosing a data repository as described above.

[https://www.re3data.org/](https://www.re3data.org/)

************************************************

                 {{6}}
************************************************

**Some well-known and more generic repositories**

************************************************

                 {{6-7}}
************************************************

* [Zenodo](https://zenodo.org/) – a repository that enables researchers, scientists, EU projects and institutions to share and showcase multidisciplinary research results (data and publications) that are not part of the existing institutional or subject-based repositories of the research communities;
* [Dryad](http://www.datadryad.org/) – a curated general-purpose repository that makes the data underlying scientific publications discoverable, freely reusable and citable. Dryad has integrated data submission for a growing list of journals;
* [Open Science Framework (OSF)](https://osf.io/) - a scholarly commons to connect the entire research cycle. It is part network of research materials, part version control system, and part collaboration software;
* [Figshare](https://figshare.com/) – a repository that allows researchers to publish all of their research outputs in an easily citable, sharable and discoverable manner.

### Explore data repositories

You have just learned about the existence of a global registry of research data repositories that covers repositories from different academic disciplines.

Re3data.org makes it possible to search for a repository that meets your criteria.

Go to [www.re3data.org/search](http://www.re3data.org/search) and find a repository that meets all three of the following criteria:

* Certificate → CoreTrustSeal
* Data licenses → CC0 (Creative Commons 0)
* Persistent identifier (PID systems) → DOI (Digital Object Identifier)

Make use of the filters offered on the left side of the screen, as visualized here:

TODO: quiz with ELIXIR resources

************************************************

### Give clarity with (Creative Commons) licenses

In order to publish your data and make it reusable, you require a license. A license creates clarity and certainty for potential users of your data. A license is not an option for all data; some of it may be too confidential or privacy-sensitive to be published.

**Creative Commons licenses**

Licenses such as the [Creative Commons](https://creativecommons.org/share-your-work/licensing-types-examples/) (CC) licenses replace 'all rights reserved' copyright with 'some rights reserved'. There are seven standard CC licenses. CC-BY is the most commonly used license, in which attribution is mandatory when using data. You can also choose restrictions like non-commercial, no derivatives, or share alike. Creative Commons offers a [guide](https://creativecommons.org/choose/?lang=en) to help you determine your preferred license.

![Creative Commons](img/CC.png)

**Assigning a license to your data**

Assigning licenses to data can also have disadvantages. Licenses are static and do not change with the quick developments in the field of research data. Therefore, some data repositories work with a CC0 license whereby no rights are reserved. Instructions regarding use are completed with codes of conduct, which may be adapted more easily.

A short movie explaining the different Creative Commons elements is shown below. Remember that sharing without a license can still lead to conflicts.

TODO: add video on CC licenses?

**Question**

We are very interested to know what license you would choose if you were to share the underlying research data of your most recent publication.  

An explanation for each license can be found by clicking on the links below.

1. CC BY: [Attribution](https://creativecommons.org/share-your-work/licensing-types-examples/licensing-examples/#by)
2. CC BY-SA: [Attribution ShareAlike](https://creativecommons.org/share-your-work/licensing-types-examples/licensing-examples/#sa)
3. CC BY-ND: [Attribution-NoDerivs](https://creativecommons.org/share-your-work/licensing-types-examples/licensing-examples/#nd)
4. CC BY-NC: [Attribution-NonCommercial](https://creativecommons.org/share-your-work/licensing-types-examples/licensing-examples/#nc)
5. CC BY-NC-SA: [Attribution-NonCommercial-ShareAlike](https://creativecommons.org/share-your-work/licensing-types-examples/licensing-examples/#by-nc-sa)
6. CC BY-NC-ND: [Attribution-NonCommercial-NoDerivs](https://creativecommons.org/share-your-work/licensing-types-examples/licensing-examples/#by-nc-nd)
7. CC0: [Public Domain](https://creativecommons.org/share-your-work/public-domain/)

### Publishing in a data journal

Data journals are publications whose primary purpose is to publish datasets. They enable you as an author to focus on the data itself, rather than producing an extensive analysis of the data which occurs in the traditional journal model. Fundamentally, data journals seek to:

* Promote scientific accreditation and reuse;
* Improve transparency of scientific methods and results;
* Support good data management practices;
* Provide an accessible and permanent route to the dataset.

**The benefits of publishing in a data journal**

Publishing in a data journal may be of interest to researchers and data producers for whom data is a primary research output. In some cases, the publication cycle may be quicker than that of traditional journals, and where there is a requirement to deposit data in an "approved repository", long-term curation and access to the data is assured.

Publishing a data paper may be regarded as best practice in data management as it:

* Includes an element of peer review of the dataset;
* Maximises opportunities for reuse of the dataset;
* Provides academic accreditation for data scientists as well as for front-line researchers.
(source: [ANDS Guide](http://www.ands.org.au/working-with-data/publishing-and-reusing-data/data-journals))

**General and disciplinary data journals**

There are data journals for various disciplines and also more general data journals exist. A widespread standard PID is the DOI. DOI stands for ‘Digital Object Identifier’. A DOI is an alphanumeric string assigned to an object which allows for an object to be identified over time. Often a DOI will be presented as a link which looks like: https://doi.org/10.1109/5.771073. There are other identifiers available which some repositories may use instead. If you are depositing in a reputable repository then you should be given some type of persistent identifier which you can use to cite and link to your data.

Examples of generic data journals:

* [Scientific Data](http://www.nature.com/sdata/about)  
* [Data in Brief](http://www.journals.elsevier.com/data-in-brief)   
* [Data Science Journal](http://www.codata.org/publications/data-science-journal)

Examples of disciplinary data journals:

TODO: check for life science additions

Open archaeology data;
Earth System Science Data;
Research Data Journal for the Humanities and Social Sciences.

### How to cite a dataset

Citations to your data can add to your academic impact.

A citation should include enough information so that the exact version of the data being cited can be located. Including a Persistent Identifier (PID) in the citation ensures that even if the location of the data changes, the PID will always link to the data that were used.

You can indicate in your (Creative Commons) license or user agreement that you want your data cited when reused.

Data citations work just like book or journal article citations and can include the following information:

* Author;
* Year;
* Dataset title;
* Repository;
* Version;
* Persistent IDentifier (PID), often works as a functional link/URL.

**Examples**

A widespread standard PID is the DOI. DOI stands for ‘Digital Object Identifier’. A DOI is an alphanumeric string assigned to an object which allows for an object to be identified over time. Often a DOI will be presented as a link which looks like: https://doi.org/10.1109/5.771073. There are other identifiers available which some repositories may use instead. If you are depositing in a reputable repository then you should be given some type of persistent identifier which you can use to cite and link to your data.

Irino, T; Tada, R (2009): Chemical and mineral compositions of sediments from ODP Site 127‐797. Geological Institute, University of Tokyo. http://dx.doi.org/10.1594/PANGAEA.726855


**Tips**

Tip1: Get a PID at the data repository of your choice.
Tip2: Is your PID a DOI and do you want to cite it in the format of a specific journal? Use the [DOI formatter](https://citation.crosscite.org/) from CrossCite.


TODO: add short quiz

### FAIR data

FAIR stands for ‘Findable, Accessible, Interoperable, and Reusable’. The FAIR data principles act as an international guideline for the result of high-quality data management.

With the increase in volume, complexity and creation speed of data, humans are more and more relying on computational support for dealing with data. The principles were defined with the focus on machine-actionability, i.e. the capacity of computational systems to find, access, interoperate and reuse data with none or minimal human intervention.

* F – Findable

By using correct metadata to describe the data, it will be findable. By using a persistent identifier the data can be found by computer systems automatically.

* A – Accessible

The data should be accessible for the long term. Even when underlying data is not accessible, the describing metadata should remain available.

* I – Interoperable

The data can be used and combined with other datasets. To achieve this, the data should be stored in generic file types, not in software specific file types.

* R – Reusable

The options for reuse should be stated clearly in a license. Without a license there is no certainty about the options for reuse and creator rights are implicit.

**How to achieve FAIR data**

In general, having a good data management plan will lead to FAIR data. In the case of privacy-sensitive data, it is possible to meet the criteria, but not to share the data openly. In this case you can make sure that a well-described dataset can be found online, while preventing the underlying data to be downloaded and used without permission.

If you anonymise your data, presuming the data is of limited sensitivity and you are very sure the data cannot lead back to the persons involved, you can share your data openly.

The FAIR Guiding Principles were put together and published in Scientific Data (Mark D. Wilkinson et al., “The FAIR Guiding Principles for Scientific Data Management and Stewardship,” Scientific Data 3 (March 15, 2016): 160018.).

TODO: add question H5P quiz?

### Open science

“Open Science is the practice of science in such a way that others can collaborate and contribute, where research data, lab notes and other research processes are freely available, under terms that enable reuse, redistribution and reproduction of the research and its underlying data and methods.”

(Source: ]FOSTER](https://www.fosteropenscience.eu/foster-taxonomy/open-science-definition)).

You have learned that good data management contributes to the findability, accessibility, interoperability and reusability of your research data. This does not necessarily mean that you should make your data openly available. But to open up data, you do need good data management from the earliest possible stage of your research project.

TODO: add links to ORION course or other relevant elements
Flemish open science plan?

### Write your data management plan for your data reuse

Go to DMPonline and open your draft data management plan created in the Introduction.

You have now completed the module on data sharing and availability for reuse. You should be able to complete the following questions in the section ‘Data availability for reuse’:

* What secondary use of your data is intended or foreseeable?
* Where will you make your data available?
* What access and usage conditions will apply?

## Rounding up

### Introduction to rounding up

You have almost reached the end of this course on research data management.

You have learned about data collection, data documentation, data storage and security, selection and preservation and making data available for reuse.

We are very curious to know if this course has helped you write your Data Management Plan (DMP).

To round up:

* We want to remind you of the DMP review service of RDM Support;
* We want to share some good practices of data management with you;
* We are collecting experiences of researchers at Utrecht University of how a DMP has helped them with their research. Please read their stories and/or share your story;
* We invite you to fill out the evaluation of this online training. This will help us to further develop this training and future learners can benefit from this. Thank you very much!

### DMP review service

You can have your data management plan (DMP) checked by the specialists of Research Data Management Support. You can get in touch if you are unsure about sections in your DMP or when you doubt whether your plan fits the requirements of your research funder.

When you are in the process of writing a proposal for a research funder and you want a check on the data section, you can also contact the Research Support Office (RSO) of your faculty.

### Researchers sharing their experiences

TODO: add stories if available or links to resources

TODO: merge with experiences

### More data stories

Challenges in irreproducible research
special issue in Nature, 7 Oct 2015

There is growing alarm about results that cannot be reproduced.  Explanations include increased levels of scrutiny, complexity of experiments and statistics, and pressures on researchers. Journals, scientists, institutions and funders all have a part in tackling reproducibility. Nature has taken substantive steps to improve the transparency and robustness in what they publish, and to promote awareness within the scientific community.

Data stories in environmental science
collected by DataONE

Success stories and cautionary tales from researchers related to their experiences with managing and sharing scientific research data as collected by DataONE.

Advantages of data sharing
by John-Alan Pascoe of Delft University of Technology

John-Alan Pascoe, researcher at the Faculty of Aerospace Engineering at Delft University of Technology, explains the advantages he experienced after sharing his raw and derived data in the data archive of 4TU.ResearchData.

!?[John-Alan Pascoe](https://youtu.be/Q7vC0v988R4)

### Evaluation of training

TODO: link to questionnaire
