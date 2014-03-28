'use strict';


// Declare app level module which depends on filters, and services
angular.module('myApp', [
    'ngRoute',
    'myApp.filters',
    'myApp.services',
    'myApp.directives',
    'myApp.controllers'
]).
        config(['$routeProvider', function($routeProvider) {
                $routeProvider.when('/introduction', {templateUrl: 'partials/introduction.html', controller: 'MyCtrl1'});
                $routeProvider.when('/taxform', {templateUrl: 'partials/taxform.html', controller: 'MyCtrl1'});
                $routeProvider.when('/datainfo', {templateUrl: 'partials/datainfo.html', controller: 'MyCtrl1'});
                $routeProvider.when('/datacreator', {templateUrl: 'partials/datacreator.html', controller: 'MyCtrl1'});
                $routeProvider.when('/geneassignment', {templateUrl: 'partials/geneassignment.html', controller: 'MyCtrl1'});
                $routeProvider.when('/genedatainfo', {templateUrl: 'partials/genedatainfo.html', controller: 'MyCtrl1'});
                $routeProvider.when('/summaryform', {templateUrl: 'partials/summaryform.xhtml', controller: 'MyCtrl1'});
                $routeProvider.when('/server_crew', {templateUrl: 'partials/server_crew.html', controller: 'MyCtrl1'});
                $routeProvider.otherwise({redirectTo: '/introduction'});
            }]);
