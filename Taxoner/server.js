/*
 * Dependencies
 */
var express = require('express'),
        app = express(),
        http = require('http'),
        url = require('url'),
        exec = require('child_process').exec,
        libpath = require('path'),
        fs = require('fs'),
        mime = require('mime'),
        readline = require('readline'),
        stream = require('stream');
/*
 * Server Config
 */
var host = "127.0.0.1",
        htmlPort = "8081",
        filePort = "8082",
        cmdPort = "8083",
        statusPort = "8084",
        cmdURL = "http://" + host + ":" + cmdPort,
        htmlURL = "http://" + host + ":" + htmlPort,
        fileURL = "http://" + host + ":" + filePort,
        statusURL = "http://" + host + ":" + statusPort;

var nodesMap = {};

var header = "<!doctype html>"
        + "<html lang='en'>"
        + "<head>"
        + "<meta charset='utf-8'>"
        + "<link rel='stylesheet' type='text/css' href='" + htmlURL + "/css/layout.css' /> "
        + "<link rel='stylesheet' type='text/css' href='" + htmlURL + "/css/style.css' />"
        + "<link rel='stylesheet' type='text/css' href='" + htmlURL + "/css/colors.css' />"
        + "</style>"
        + "</head>"
        + "<body>";

var end = "</body></html>";
/*
 * Server to run commands over http
 */
http.createServer(function(req, res) {
    var parsedUrl = url.parse(req.url, true);
    var cmd = parsedUrl.query['cmd'];
    var async = parsedUrl.query['async'];
    res.writeHead(200, {'Content-Type': 'text/plain'});
    if (cmd) {
        var child = exec(cmd, function(error, stdout, stderr) {
            var result = '{"stdout":' + stdout + ',"stderr":"' + stderr + '","cmd":"' + cmd + '"}';
            res.end(stdout);
        });
    } else {
        var result = '{"stdout":"' + '' + '","stderr":"' + 'cmd is mandatory' + '","cmd":"' + cmd + '"}';
        res.end(result + '\n');
    }
    if (async === "true") {
        var result = '{"stdout":"async request' + '' + '","stderr":"' + '' + '","cmd":"' + cmd + '"}';
        res.end(result + '\n');
    }

}).listen(cmdPort, host);
console.log('CMD server running at ' + cmdURL);

/*
 * Server to print the Log file
 */
http.createServer(function(req, res) {
    var parsedUrl = url.parse(req.url, true);
    var filename = parsedUrl.query['file'];
    var n = parsedUrl.query['n'];
    var h;

    if (!n || n === 2) {
        n = 2;
        h = 1;
    } else {
        h = n;
    }

    if (filename) {
        fs.exists(filename, function(exists) {
            if (!exists) {
                res.writeHead(404, {
                    "Content-Type": "text/plain"
                });
                res.write("404 Not Found\n");
                res.end();
                return;
            }

            if (fs.statSync(filename).isFile) {
                res.writeHead(200, {
                    "Content-Type": "text/html"
                });
                var cmd = "cat " + filename;
                if (n != 0) {
                    cmd = cmd + " | sed 's/\\r/\\n/g' | tail -n " + n + " | head -n " + h;
                }
                var child = exec(cmd, function(error, stdout, stderr) {
                    var result = header
                            + "<script>var a = window.setTimeout('window.location.reload();', 2000);</script>"
                            + "<div align='center'>"
                            + "<h1 class='title' id='page-title'>Log viewer</h1>"
                            + "<a href='" + cmdURL + "/?cmd=xdg-open%20" + filename.replace('stdout.log', '') + "' target='_blank'/>Click to go to your output dir</a><br>"
                            + "</div>"
                            + "<br><br>";

                    if (typeof (stdout) === "undefined" || stdout === null || stdout === "") {
                        result = result + "Wait .... ";
                    } else {
                        result = result + "<pre>" + stdout + "</pre>";
                    }

                    result = result
                            + end;

                    res.end(result);
                });
            }
        });
    } else {
        var result = "This is a server to list files\n" +
                "URL: "
                + fileURL + "/?file=filename.txt";
        res.end(result + '\n');
    }
}).listen(filePort, host);
console.log('File server running at ' + fileURL);

/*
 * Server to show status
 */
http.createServer(function(req, res) {
    var parsedUrl = url.parse(req.url, true);
    var filename = parsedUrl.query['file'];
    var tax = parsedUrl.query['tax'];
    res.writeHead(200, {'Content-Type': 'text/html'});

    setNodesRank(setNodesNames, filename, tax, res);

}).listen(statusPort, host);
console.log('Summary server running at ' + statusURL);

/*
 * Server for the http site
 */
app.configure(function() {
    app.use(express.static(__dirname + '/app'));
    app.use(express.logger('dev'));
    app.use(express.bodyParser());
    app.use(express.methodOverride());
});
app.listen(htmlPort);
console.log('HTML server running at ' + htmlURL);

/**
 * This function create the server who parse the results files for Taxoner and 
 * GeneAssignment to create the summary tables
 * 
 * @param {type} filename the file to be parsed
 * @param {type} res the html response
 */
function server(filename, tax, res) {
    if (typeof String.prototype.endsWith !== 'function') {
        String.prototype.endsWith = function(suffix) {
            return this.indexOf(suffix, this.length - suffix.length) !== -1;
        };
    }
    if (filename) {
        fs.exists(filename, function(exists) {
            if (!exists) {
                var result = header
                        + "<script>var a = window.setTimeout('window.location.reload();', 2000);</script>"
                        + "<h1 class='title' id='page-title'>Summary</h1>"
                        + "Waiting for the result file";
                +end;
                res.end(result);
                return;
            } else {
                var file = fs.statSync(filename);
                if (file.isFile) {
                    if (file["size"] === 0) {
                        var result = header
                                + "<script>var a = window.setTimeout('window.location.reload();', 2000);</script>"
                                + "<h1 class='title' id='page-title'>Summary</h1>"
                                + "Waiting for the result file";
                        +end;
                        res.end(result);
                        return;
                    } else {
                        if (tax === 'true') {
                            var instream = fs.createReadStream(filename);
                            var outstream = new stream;
                            outstream.readable = true;
                            outstream.writable = true;

                            var result = header;
                            var taxs = {};

                            var rl = readline.createInterface({
                                input: instream,
                                output: outstream,
                                terminal: false
                            });

                            rl.on('line', function(line) {
                                fields = line.split("\t");
                                if (fields.length > 1) {
                                    if (typeof (taxs[fields[1]]) === "undefined") {
                                        taxs[fields[1]] = {};
                                        taxs[fields[1]]['taxId'] = fields[1];
                                        taxs[fields[1]]['total'] = 1;
                                        taxs[fields[1]]['gi'] = {};
                                        taxs[fields[1]]['gi'][fields[2]] = 1;
                                    } else {
                                        if (typeof (taxs[fields[1]]['gi'][fields[2]]) === "undefined") {
                                            taxs[fields[1]]['gi'][fields[2]] = 1;
                                        } else {
                                            taxs[fields[1]]['gi'][fields[2]] = taxs[fields[1]]['gi'][fields[2]] + 1;
                                        }
                                        taxs[fields[1]]['total'] = taxs[fields[1]]['total'] + 1;
                                    }
                                }
                            }).on('close', function() {
                                var sortable = [];
                                for (var key in taxs) {
                                    sortable.push([key, taxs[key]['total']]);
                                }
                                sortable.sort(function(a, b) {
                                    return b[1] - a[1];
                                });
                                result = result
                                        + "<h1 class='title' id='page-title'>Summary</h1>"
                                        + "<table>"
                                        + "<tr><th>Taxonomy</th>"
                                        + "<th>Rank</th>"
                                        + "<th>No. of Reads</th></tr>";
                                for (var i = 0; i < sortable.length; i++) {
                                    var key = sortable[i][0];
                                    var taxId = taxs[key]['taxId'];
                                    if (typeof (nodesMap[taxId]) === "undefined") {
                                        result = result
                                                + "<tr><td>"
                                                + "<a href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="
                                                + taxId
                                                + "' target='_blank'>"
                                                + " ("
                                                + taxId
                                                + ")"
                                                + "</a>"
                                                + "</td><td>"
                                                + "</td><td>"
                                                + taxs[taxId]['total']
                                                + "</td></tr>";
                                    } else {
                                        result = result
                                                + "<tr><td>"
                                                + "<a href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="
                                                + taxId
                                                + "' target='_blank'>"
                                                + nodesMap[taxId]['name']
                                                + " ("
                                                + taxId
                                                + ")"
                                                + "</a>"
                                                + "</td><td>"
                                                + nodesMap[taxId]['rank']
                                                + "</td><td>"
                                                + taxs[taxId]['total']
                                                + "</td></tr>";
                                    }
                                }
                                result = result
                                        + "</table>";
                                res.end(result);
                            });
                        } else {
                            var instream = fs.createReadStream(filename);
                            var outstream = new stream;
                            outstream.readable = true;
                            outstream.writable = true;

                            var result = header;
                            var taxs = {};

                            var rl = readline.createInterface({
                                input: instream,
                                output: outstream,
                                terminal: false
                            });

                            rl.on('line', function(line) {
                                fields = line.split("\t");
                                if (fields.length > 1) {
                                    if (typeof (taxs[fields[1]]) === "undefined") {
                                        taxs[fields[1]] = {};
                                        taxs[fields[1]]['taxId'] = fields[1];
                                        taxs[fields[1]]['total'] = parseInt(fields[4]);
                                        taxs[fields[1]]['prot'] = 1;
                                        taxs[fields[1]]['gi'] = {};
                                        taxs[fields[1]]['gi'][fields[2]] = parseInt(fields[4]);
                                    } else {
                                        if (typeof (taxs[fields[1]]['gi'][fields[2]]) === "undefined") {
                                            taxs[fields[1]]['gi'][fields[2]] = parseInt(fields[4]);
                                        } else {
                                            taxs[fields[1]]['gi'][fields[2]] = taxs[fields[1]]['gi'][fields[2]] + parseInt(fields[4]);
                                        }
                                        taxs[fields[1]]['total'] = taxs[fields[1]]['total'] + parseInt(fields[4]);
                                        taxs[fields[1]]['prot'] = taxs[fields[1]]['prot'] + 1;
                                    }
                                }
                            }).on('close', function() {
                                var sortable = [];
                                for (var key in taxs) {
                                    sortable.push([key, taxs[key]['total']]);
                                }
                                sortable.sort(function(a, b) {
                                    return b[1] - a[1];
                                });
                                result = result
                                        + "<h1 class='title' id='page-title'>Summary</h1>"
                                        + "<table>"
                                        + "<tr><th>Taxonomy</th>"
                                        + "<th>Rank</th>"
                                        + "<th>No. of Genes</th>"
                                        + "<th>No. of Reads</th></tr>";
                                for (var i = 0; i < sortable.length; i++) {
                                    var key = sortable[i][0];
                                    var taxId = taxs[key]['taxId'];
                                    if (typeof (nodesMap[taxId]) === "undefined") {
                                        console.log("TaxId: ", taxId);
                                        result = result
                                                + "<tr><td>"
                                                + "<a href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="
                                                + taxId
                                                + "' target='_blank'>"
                                                + " ("
                                                + taxId
                                                + ")"
                                                + "</a>"
                                                + "</td><td>"
                                                + "</td><td>"
                                                + taxs[taxId]['prot']
                                                + "</td><td>"
                                                + taxs[taxId]['total']
                                                + "</td></tr>";
                                    } else {
                                        result = result
                                                + "<tr><td>"
                                                + "<a href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="
                                                + taxId
                                                + "' target='_blank'>"
                                                + nodesMap[taxId]['name']
                                                + " ("
                                                + taxId
                                                + ")"
                                                + "</a>"
                                                + "</td><td>"
                                                + nodesMap[taxId]['rank']
                                                + "</td><td>"
                                                + taxs[taxId]['prot']
                                                + "</td><td>"
                                                + taxs[taxId]['total']
                                                + "</td></tr>";
                                    }
                                }
                                result = result
                                        + "</table>";
                                res.end(result);
                            });
                        }
                    }
                }

            }
        });
    } else {
        var result = header
                + "<script>var a = window.setTimeout('window.location.reload();', 2000);</script>"
                + "This server analyze the result files<br>"
                + "The useage:<br>"
                + statusURL + "/?file=filename.txt"
                + end;
        res.end(result);
    }
}

/**
 * This function parse the nodes.dmp file and creates a nodesMap object to 
 * storage the data
 * 
 * @param {type} callback The clallback function to execute after parse the file
 * @param {type} filename the result file to be parsed
 * @param {type} res 
 * @returns {undefined}
 */
function setNodesRank(callback, filename, tax, res) {
    if (typeof (nodesMap[2]) === "undefined") {
        console.log("Reading the nodes file from databases/nodes.dmp");
        var instream = fs.createReadStream("databases/nodes.dmp");
        var outstream = new stream;
        outstream.readable = true;
        outstream.writable = true;

        var rl = readline.createInterface({
            input: instream,
            output: outstream,
            terminal: false
        });

        rl.on('line', function(line) {
            fields = line.split("|");
            if (fields.length > 1) {
                nodesMap[fields[0].trim()] = {};
                nodesMap[fields[0].trim()]['rank'] = fields[2].trim();
            }
        }).on('close', function() {
            callback(server, filename, tax, res);
        });
    } else {
        callback(server, filename, tax, res);
    }
}

function setNodesNames(callback, filename, tax, res) {
    if (typeof (nodesMap[2]['name']) === "undefined") {
        console.log("Reading the name file from databases/names");
        var instream = fs.createReadStream("databases/names");
        var outstream = new stream;
        outstream.readable = true;
        outstream.writable = true;

        var rl = readline.createInterface({
            input: instream,
            output: outstream,
            terminal: false
        });

        rl.on('line', function(line) {
            fields = line.split("\t");
            if (fields.length > 1) {
                if (typeof (nodesMap[fields[0].trim()]['rank']) !== "undefined") {
                    nodesMap[fields[0].trim()]['name'] = fields[1].trim();
                }
            }
        }).on('close', function() {
            callback(filename, tax, res);
        });
    } else {
        callback(filename, tax, res);
    }
}

