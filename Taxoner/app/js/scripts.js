/**
 * Function to upload the html page to start the command
 * 
 * @param {String} strURL the URL to be load
 * @param {String} formName the Form name
 * @param {String} logLines number of lines in the log
 */
function xmlhttpUpload(strURL, formName, logLines) {
    var xmlreq = new XMLHttpRequest();
    var result = getFormData(strURL, formName);
    if (result == null) {
        return null;
    }
    console.log('url:', result[0]);
    console.log('outdir:', result[1]);
    xmlreq.open('POST', result[0], true);
    xmlreq.onreadystatechange = function() {
        if (xmlreq.readyState == 4) {
            var html =
                    '<iframe id="sumFrame" width="100%" seamless="seamless" '
                    + 'height="340px" scrolling="yes" frameborder="0" src="http://127.0.0.1:8084/?file='
                    + result[2] + result[3]
                    + '"></iframe><br>'
                    + '<iframe id="extFrame" width="100%" seamless="seamless" '
                    + 'height="340px" scrolling="yes" frameborder="0" src="http://127.0.0.1:8082/?file='
                    + result[1] + '.log&n=' + logLines
                    + '"></iframe><br>';
            updateElementById("result", html);
        }
    };
    xmlreq.send(result[0]);
}

/**
 * Function to upload the html page to start the command
 * 
 * @param {type} strURL the URL to be load
 * @param {type} formName the Form name
 */
function xmlhttpUpload3(strURL, formName) {
    var xmlreq = new XMLHttpRequest();
    var result = getFormData(strURL, formName);
    if (result == null) {
        return null;
    }
    console.log('url:', result[0]);
    xmlreq.open('POST', result[0], true);
    xmlreq.onreadystatechange = function() {
        if (xmlreq.readyState === 4) {
            var html =
                    '<iframe id="sumFrame" width="100%" seamless="seamless" '
                    + 'height="1000px" scrolling="yes" frameborder="0" src="'
                    + result[0]
                    + '"></iframe;';
            updateElementById("result", html);
        }
    };
    xmlreq.send(result[0]);
}

/**
 * Extrat the form data and form the POST URL
 * 
 * @param {String} strURL the basic URL to be add the data
 * @param {String} formName the Form name
 * @returns {Array} array of 2 elements. [0] the full URL to be upload, 
 * [1] the command output file
 */
function getFormData(strURL, formName) {
    var result = new Array();
    if (formName == "taxForm") {
        var seq = document.getElementById("seq").value;
        var fasta = document.getElementById("fasta").checked;
        var megan = document.getElementById("megan").checked;
        var o = document.getElementById("o").value;
        var dbPath = document.getElementById("dbPath");
        var taxpath = document.getElementById("taxpath").value;
        var dbList = document.getElementById("dbList").value;
        var dbNames = document.getElementById("dbNames");
        var nohostfilter = document.getElementById("nohostfilter").checked;
        var bt2maxhits = document.getElementById("bt2maxhits").value;
        var neighborscore = document.getElementById("neighborscore").value;
        var bt2allhits = document.getElementById("bt2allhits").checked;
        var onlyneighbor = document.getElementById("onlyneighbor").checked;
        var p = document.getElementById("p").value;
        var paired = document.getElementById("paired").value;
        var I = document.getElementById("I").value;
        var X = document.getElementById("X").value;

        if (typeof (seq) == "undefined" || seq == null || seq == "") {
            alert("Set the input file (reads), please");
            return null;
        }

        if (typeof (o) == "undefined" || o == null || o == "") {
            alert("Set the output folder path, please");
            return null;
        }

        if ((typeof (dbList) == "undefined" || dbList == null || dbList == "") &&
                (typeof (dbPath) == "undefined" || dbPath == null || dbPath == "")) {
            alert("Set the path to the bowtie2 indexes of database, please");
            return null;
        }

        if (typeof (taxpath) == "undefined" || taxpath == null || taxpath == "") {
            alert("Set the path to nodes.dmp file from NCBI Taxonomy database, please");
            return null;
        }

        result[0] = strURL + "%20-v%20%20-seq%20" + seq
                + "%20-o%20" + o
                + "%20-taxpath%20" + taxpath;

        if (fasta == true) {
            result[0] = result[0] + "%20-fasta";
        }

        if (megan == true) {
            result[0] = result[0] + "%20-megan";
        }

        if (dbList == 0) {
            /* Align with all DB */
            result[0] = result[0] + "%20-dbPath%20databases/bowtie2";
        } else if (dbList == 1) {
            /* Align with Bacteria */
            result[0] = result[0] + "%20-dbPath%20databases/bowtie2"
                    + "%20-dbNames%20bact.0.fasta,bact.1.fasta,bact.2.fasta,bact.3.fasta";
            //+ "%20-dbNames%20bact.0.fasta";
        } else if (dbList == 2) {
            /* Align with Archaea */
            result[0] = result[0] + "%20-dbPath%20databases/bowtie2"
                    + "%20-dbNames%20archaea.0.fasta";
        } else if (dbList == 3) {
            /* Align with Viruses */
            result[0] = result[0] + "%20-dbPath%20databases/bowtie2"
                    + "%20-dbNames%20virus.0.fasta";
        } else if (dbList == 4) {
            /* Align with Fungi */
            result[0] = result[0] + "%20-dbPath%20databases/bowtie2"
                    + "%20-dbNames%20fungi.0.fasta";
        } else if (dbList == 5) {
            /* Align with own DB */
            result[0] = result[0] + "%20-dbPath%20" + dbPath.value;
        }

        if (typeof (dbNames) != "undefined" && dbNames != null && dbNames != "") {
            result[0] = result[0] + "%20-dbNames%20" + dbNames.value;
        }

        if (nohostfilter == true) {
            result[0] = result[0] + "%20-no-host-filter";
        }

        if (typeof (bt2maxhits) != "undefined" && bt2maxhits != null && bt2maxhits != "") {
            result[0] = result[0] + "%20-bt2-maxhits%20" + bt2maxhits;
        }

        if (bt2allhits == true) {
            result[0] = result[0] + "%20-bt2-allhits";
        }

        if (typeof (neighborscore) != "undefined" && neighborscore != null && neighborscore != "") {
            result[0] = result[0] + "%20-neighbor-score%20" + neighborscore;
        }

        if (bt2allhits == true) {
            result[0] = result[0] + "%20-bt2-allhits";
        }

        if (onlyneighbor == true) {
            result[0] = result[0] + "%20-only-neighbor ";
        }

        if (typeof (p) != "undefined" && p != null && p != "") {
            result[0] = result[0] + "%20-p%20" + p;
        }

        if (typeof (paired) != "undefined" && paired != null && paired != "") {
            result[0] = result[0] + "%20-paired%20" + paired;
        }

        if (typeof (I) != "undefined" && I != null && I != "") {
            result[0] = result[0] + "%20-I%20" + I;
        }

        if (typeof (X) != "undefined" && X != null && X != "") {
            result[0] = result[0] + "%20-X%20" + X;
        }

        result[0] = result[0]
                + "%20>%20" + o + "/stdout.log"
                + "&async=true";

        result[1] = o + "/stdout";
        result[2] = o + "/Results/Taxonomy.txt";
        result[3] = "&tax=true";
    } else if (formName == "geneForm") {
        var inFile = document.getElementById("inFile").value + 'Results/Taxonomy.txt';
        var outDir = document.getElementById("inFile").value + 'Results';
        var binDB = document.getElementById("binDB").value;
        var indDB = document.getElementById("indDB").value;

        if (typeof (inFile) == "undefined" || inFile == null || inFile == "") {
            alert("Set the input file, please");
            return null;
        }

        if (typeof (outDir) == "undefined" || outDir == null || outDir == "") {
            alert("Set the output file, please");
            return null;
        }

        if (typeof (binDB) == "undefined" || binDB == null || binDB == "") {
            alert("Set the database binary file folder path, please");
            return null;
        }

        if (typeof (indDB) == "undefined" || indDB == null || indDB == "") {
            alert("Set the database index file, please");
            return null;
        }

        result[0] = strURL + "%20-v%20%20-t%20" + inFile
                + "%20-o%20" + outDir + "/genes.txt"
                + "%20-b%20" + binDB
                + "%20-i%20" + indDB
                + ">" + outDir + "/genes.log"
                + "&async=true";
        result[1] = outDir + "/genes";
        result[2] = outDir + "/genes.txt";
        result[3] = "&tax=false";
    } else if (formName == "sumForm") {
        var filename = document.getElementById("filename").value;
        var tax = document.getElementById("tax").checked;

        result[0] = strURL + "/?file=" +filename;

        if (tax == true) {
            result[0] = result[0] + "&tax=true";
        } else {
            result[0] = result[0] + "&tax=false";
        }
    }
    return result;
}

/**
 * Update a html element by Id
 * 
 * @param {String} id the id of the html element
 * @param {String} str the HTML code to be inserted
 */
function updateElementById(id, str) {
    document.getElementById(id).innerHTML = str;
}

/**
 * Reload the form to reset the values
 * 
 * @param {String} newhash the path to be reaload
 */
function reloadForm(newhash) {
    window.location.hash = "#" + newhash;
    window.location.reload(true);
}

/**
 * Duplicate to filds in a Form
 * 
 * @param {Form} formName the Form with the fields to be duplicated
 */
function duplicateFormInput(formName) {
    if (formName.name == "geneForm") {
        formName.outDir.value = formName.inFile.value + ".out";
    }
}

/**
 * Insert the input to use an extra bowtie database
 * 
 * @param {type} formName the form name
 */
function verifyDBList(formName) {
    if (formName.name == "taxForm") {
        if (formName.dbList.value == 5) {
            var htmlToInsert = "Path to bowtie2 indexes of database:<br><input id='dbPath' size='54'/><br>"
                    + "<font title='Example: 0.fasta,1.fasta of the bowtie indexes'>"
                    + "Names of index files to which reads will be aligned (OPTIONAL):</font><br>"
                    + "<input id='dbNames' size='54' title='Example: 0.fasta,1.fasta'/>";

            document.getElementById("otherDB").innerHTML = htmlToInsert;
        } else {
            document.getElementById("otherDB").innerHTML = "";
        }
    }
}

/**
 * Check the checkbox on the summary form
 * 
 * @param {type} id the element id
 */
function verifySummaryCheck(id) {
    if (id === "tax") {
        if (document.getElementById("tax").checked) {
            document.getElementById("gen").checked = false;
        } else {
            document.getElementById("gen").checked = true;
        }
    } else if (id === "gen") {
        if (document.getElementById("gen").checked) {
            document.getElementById("tax").checked = false;
        } else {
            document.getElementById("tax").checked = true;
        }
    }
}
