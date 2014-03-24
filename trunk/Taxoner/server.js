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
        mime = require('mime');
/*
 * Server Config
 */
var host = "127.0.0.1",
        htmlPort = "8081",
        filePort = "8082",
        cmdPort = "8083",
        cmdURL = "http://" + host + ":" + cmdPort,
        htmlURL = "http://" + host + ":" + htmlPort,
        fileURL = "http://" + host + ":" + filePort;
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
    if (async == "true") {
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

    if (!n || n == 2) {
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
                    var result = "<!doctype html>"
                            + "<html lang='en'>"
                            + "<head>"
                            + "<meta charset='utf-8'>"
                            + "<link rel='stylesheet' type='text/css' href='" + htmlURL + "/css/layout.css' /> "
                            + "<link rel='stylesheet' type='text/css' href='" + htmlURL + "/css/style.css' />"
                            + "<link rel='stylesheet' type='text/css' href='" + htmlURL + "/css/colors.css' />"
                            + "</head>"
                            + "<body><script>var a = window.setTimeout('window.location.reload();', 2000);</script>"
                            + "<div align='center'>"
                            + "<h1 class='title' id='page-title'>Log viewer</h1>"
                            + "<a href='" + cmdURL + "/?cmd=xdg-open%20" + filename.replace('stdout.log', '') + "' target='_blank'/>Click to go to your output dir</a>"
                            + "</div>"
                            + "<br><br>";

                    if (typeof (stdout) == "undefined" || stdout == null || stdout == "") {
                        result = result + "Wait .... ";
                    } else {
                        result = result + "<pre>" + stdout + "</pre>";
                    }

                    result = result
                            + "</body>"
                            + "</html>";

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
