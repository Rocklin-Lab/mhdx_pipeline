import glob
import os
from collections import OrderedDict

cwd=os.getcwd()

#class types are 'pName', 'rtName', 'dName', uses OrderedDict for dic
def write_nested_list(dic, class_type):
    out = ""
    count = 1
    for key in dic.keys():
        if class_type == 'rtNames':
            print_key = str(key)+"'"
        else:
            print_key = str(count)+'''&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'''+key
            
        out += '''                    <li>

                                <a class="protName '''+class_type+'''">'''+str(print_key)+'''</a>

                            <ul class="nolist">'''

        for path, species in dic[key]:
            out += '''                            <li>

                                    <a class="speciesName" onClick="(function(){document.getElementById('myFrame').src='''+"'"+path+"'"+''';return false;})();return false;">'''+species+'''</a>

                                </li>'''
        out += '''
        </ul>
        
        </li>
        '''
        count += 1
    return out

#prepare sets of protein names, species names, rt bins, and plot paths

overview_plot_path = snakemake.input[0]
all_files = snakemake.input[1:]
stripped = [file.split('/')[-1] for file in all_files]
names = sorted(set([('_').join(file.split('_')[:3]).split('.')[0] for file in stripped]))
display_names = [str(i+1)+'''&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'''+names[i] for i in range(len(names))]
rt_bins = sorted(set([int(float(file.split('_')[3])) for file in stripped]))
species_names = sorted(set([("_").join(file.split("_")[:-2]) for file in stripped]))
species_plot_paths = ["./plots/ic_time_series/html/"+species+"_time_series.html" for species in species_names]

#make dicts of: rt minute to species list, name to charge species

name_to_species = OrderedDict.fromkeys(names)
for name in name_to_species.keys():
    name_to_species[name] = []
    sub_specs = [species for species in species_names if name in species]
    for spec in sub_specs:
        path = [path for path in species_plot_paths if spec in path][0]
        name_to_species[name].append((path, spec))

rt_to_species = OrderedDict.fromkeys(rt_bins)
for rt_bin in rt_to_species.keys():
    rt_to_species[rt_bin] = []
    try:
        sub_specs = sorted([species for species in species_names if str(rt_bin) in species.split('_')[3][:len(str(rt_bin))]], key = lambda x: float(x.split("_")[-1]))
        for spec in sub_specs:
            path = [path for path in species_plot_paths if spec in path][0]
            rt_to_species[rt_bin].append((path, spec))
    except:
        import ipdb
        ipdb.set_trace()

up_to_pNames = '''
<!-- HTML -->

<html>

<body>

<meta name="viewport" content="width=device-width, initial-scale=1">

<link rel="stylesheet" href="https://fonts.googleapis.com/icon?family=Material+Icons">

<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0-alpha.6/css/bootstrap.min.css">

<script src="https://code.jquery.com/jquery-3.2.1.min.js"></script>

<script src="https://cdnjs.cloudflare.com/ajax/libs/tether/1.4.0/js/tether.min.js"></script>

<script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0-alpha.6/js/bootstrap.min.js"></script>

<!--

Layout of Sidebar nested-list:

<aside> - arbitrary tag distinct from <main> 

    <nav> - container for links
    
        <ul> list container of mainHeaders
      
            <li> mainHeader list element
      
                <a> class='mainHeader'
                |      
                <ul> 1st nested list container of protNames
      
                    <li> 1st nested list element - protNames
      
                        <a> class='protName'
                        |
                        <ul> 2nd nested list container
      
                            <li> 2nd nested list element - speciesName
      
                                <a> class='speciesName'

                            </li>
                        
                        </ul>
                    
                    </li>
                
                </ul>
            
            </li>
        
        </ul>  

-->

<main id="main">
    <iframe id="myFrame" src="'''+overview_plot_path+'''" height=100% width=100%></iframe>
</main>

<aside class="side-menu">

    <nav id="mySidebar" class="left-nav" onmouseover="toggleSidebar()" onmouseout="toggleSidebar()">

        <ul>

            <li>

                <a class="mainHeader" id="ov" onClick="(function(){document.getElementById('myFrame').src='''+"'"+overview_plot_path+"'"+''';return false;})();return false;"><i class="material-icons">info</i><span class="icon-text">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Overview</span></a>
            
            </li>

            <li>
                
                <a class="mainHeader"><i class="material-icons">assessment</i><span class="icon-text">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Plots</span></a>

                <input type='text' id='nameSearch' class='searchBar' onKeyUp="search(this)" placeholder="Name Search">
                
                <ul class="nolist">
'''

up_to_rtNames = '''                </ul>
            
            </li>

            <li>
                
                <a class="mainHeader"><i class="material-icons">query_builder</i><span class="icon-text">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Chromatogram</span></a>
                
                <ul class="nolist">'''

#                <input type='text' id='rtSearch' class='searchBar' onKeyUp="search(this)" placeholder="RT Search"> #

up_to_dNames = '''                </ul>
            
            </li>

            <li>
                
                <a class="mainHeader"><i class="material-icons">list</i><span class="icon-text">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Datasheets</span></a>

                <input type='text' id='dataSearch' class='searchBar' onKeyUp="search(this)" placeholder="Name Search">
                
                <ul class="nolist">'''

after_dNames = '''                </ul>
            
            </li>

        </ul>

    </nav>

</aside>

<!-- /HTML -->



<!-- CSS -->

<style>

a {
  text-decoration: none !important;
  font-family: Helvetica;
  
}

aside nav.left-nav {
    height: 100%;
    width: 60px;
    position: fixed;
    z-index: 1;
    top: 0;
    left: 0;
    background-color: #111;
    overflow-x: hidden;
    transition: 0.5s;
    padding-top: 0px;
    white-space: nowrap;
}

aside nav.left-nav ul {
    margin-left: 0;
    padding-left: 0
}

aside nav.left-nav ul li {
    border-bottom: 2px solid #fff;
    font-size: 1.75rem
}

aside nav.left-nav ul li a.mainHeader {
    color: #fff;
    font-weight: 500;
    display: block;
    margin-top: 0px;
    padding: 1.6rem; 
    padding-left: 1.25rem;
    font-size: 1.25rem;
    border-bottom: 3px solid #fff;
}

aside nav.left-nav ul li a.mainHeader:hover {
    color: #fff;
    background-color: #17629e;
    text-decoration: none
}

aside nav.left-nav ul li ul {
    display: none;
}

.searchBar {
    height: 40px;
    width: 190px;
    margin-left: 3.75rem;
    font-size: 1.25rem;
    display: none;
}

aside nav.left-nav ul li ul li a.protName {
    font-size: 1.25rem;
    font-weight: 200;
    padding-top: 0.6rem;
    padding-bottom: 0.55rem;
    padding-left: 1rem;
    padding-right: 14rem;
    color: #fff;
}

aside nav.left-nav ul li ul li a.protName.active {
    background-color: #fff;
    color: #17629e;
    
}

aside nav.left-nav ul li ul li a.protName:hover {
    color: #fff;
    background-color: #17629e

}

aside nav.left-nav ul li ul li ul{
    display: none;
}

aside nav.left-nav ul li ul li ul li a.speciesName {
    font-size: 0.75rem;
    font-weight: 100;
    padding-top: 0.85rem;
    padding-bottom: 0.5rem;
    padding-left: 4rem;
    padding-right: 10rem;
    color: #fff;
}

aside nav.left-nav ul li ul li ul li a.speciesName.active {
    background-color: #fff;
    color: #17629e;
}

aside nav.left-nav ul li ul li ul li a.speciesName:hover {
    color: #fff;
    background-color: #17629e
}

#main {
    transition: margin-left .5s;
    margin-left: 60px;
}

.material-icons,
    .icon-text {
        vertical-align: middle;
    }
        
.material-icons {
    padding-bottom: 0px;
}

</style>

<!-- /CSS -->



<!-- JavaScript -->

<script>

    /* Slides sidebar mainHeader down and sets to active on 1st click */
    $('.side-menu nav ul li').on('click', 'a.mainHeader', function(e){
        if ($(this).parent().children('ul').length){
            e.preventDefault();
            $(this).addClass('active');
            $(this).parent().children('ul').slideDown();
            $(this).parent().children('input.searchBar').slideDown();
        }       
    });

    /* Slides sidebar mainHeader up and sets to inactive on 2nd click */
    $('.side-menu nav ul li').on('click', 'a.mainHeader.active', function(e){
        e.preventDefault();
        $(this).removeClass('active');
        $(this).parent().children('ul').slideUp();
        $(this).parent().children('input.searchBar').slideUp();
    });

    /* Slides sidebar protNames down and sets to active on 1st click, deactivates any other protName element that doesn't have an active speciesName */
    $('.side-menu nav ul li ul li').on('click', 'a.protName', function(e){
        e.preventDefault();
        
        //Add active and slide down if len child list > 0
        $(this).addClass('active');
        if ($(this).parent().children('ul').length){
            $(this).parent().children('ul').slideDown();
        }

        

        //Look for other active protNames, if they don't have active speciesNames, slide them up and deactivate
        var outer = this; //to pass into filter function
        var activeNames = $('.protName.active').filter(function(){
            return $(this).parent().children('ul').children('li').children('.active').length == 0 && !$(this).get(0).isSameNode(outer)
        });

        if (activeNames.length > 0){
            activeNames.parent().children('ul').slideUp();
            activeNames.parent().children('.searchBar').slideUp();
            activeNames.removeClass('active');
        }
    });

    /* Slides sidebar protNames up and sets to inactive on 2nd click, or stops open header from sliding up if  */
    $('.side-menu nav ul li ul li').on('click', 'a.protName.active', function(e){
        e.preventDefault();
        if ($(this).parent().children('ul').children('li').children('.active').length > 0){
            $(this).parent().children('ul').slideDown();
        } 
        else{
            $(this).parent().children('ul').slideUp();  
            if ($(this).parent().children('ul').children('li').children('.active').length == 0){
                $(this).removeClass('active');
            }
        }
    });

    /* Sets sidebar speciesName to active on 1st click, deactivates any other active speciesName element */
    $('.side-menu nav ul li ul li ul li').on('click', 'a.speciesName', function(e){
        e.preventDefault();     
        $('.speciesName.active').removeClass('active');
        $(this).addClass('active');
    });
    
    /* Sets sidebar speciesName to inactive on 2nd click */
    $('.side-menu nav ul li ul li ul li').on('click', 'a.speciesName.active', function(e){
        e.preventDefault();
        $(this).removeClass('active');
    });

    /* On clicking Overview, deactivates all active elements, closes all expanded menus */
    $("#ov").on('click', function(e){
        e.preventDefault();

        $('.mainHeader.active').parent().children('ul').slideUp();
        $('.mainHeader.active').parent().children('input.searchBar').slideUp();
        $('.mainHeader.active').removeClass('active');

        $('.protName.active').parent().children('ul').slideUp();
        $('.protName.active').removeClass('active');

        $('.speciesName.active').removeClass('active');
    });

    /* Shows/hides sidebar on mouseover/mouseout */
    var mini = true;
    function toggleSidebar() {
        if (mini) {
            document.getElementById("mySidebar").style.width = "250px";
            document.getElementById("main").style.marginLeft = "250px";
            this.mini = false;
        } 

        else {
            document.getElementById("mySidebar").style.width = "60px";
            document.getElementById("main").style.marginLeft = "60px";
            this.mini = true;
        }
    }

    function search(ele){

        var filter, names, i, txtValue;

        filter = document.getElementById(ele.id).value.toUpperCase();

        if (ele.id == "nameSearch"){
            names = document.getElementsByClassName("pNames");
        } else{
            if(ele.id == "rtSearch"){
                names = document.getElementsByClassName("rtNames");
            } else{
                if(ele.id == "dataSearch"){
                    names = document.getElementsByClassName("dNames");
                }
            }
        }

        for (i = 0; i < names.length; i++){
            
            txtValue = names[i].innerText.toUpperCase();
            
            if (txtValue.indexOf(filter) == -1){
                names[i].parentElement.style.display = "none";
            } 
           
            else{
                names[i].parentElement.style.display = "";
            }
        }
    }

</script>

<!-- /JavaScript -->

</body>
</html>'''

#make html file

proc = ""
proc += up_to_pNames
proc += write_nested_list(name_to_species, 'pNames')
proc += up_to_rtNames
proc += write_nested_list(rt_to_species, 'rtNames')
proc += after_dNames

#save to output

with open(snakemake.output[0], 'wt') as file:
    file.write(proc)