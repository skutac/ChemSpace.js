
/*DEEP MERGE*/
function isMergeableObject(e){return e&&"object"==typeof e&&"[object RegExp]"!==Object.prototype.toString.call(e)&&"[object Date]"!==Object.prototype.toString.call(e)}function emptyTarget(e){return Array.isArray(e)?[]:{}}function cloneIfNecessary(e,r){return r&&!0===r.clone&&isMergeableObject(e)?deepmerge(emptyTarget(e),e,r):e}function defaultArrayMerge(e,r,t){var a=e.slice();return r.forEach(function(r,c){void 0===a[c]?a[c]=cloneIfNecessary(r,t):isMergeableObject(r)?a[c]=deepmerge(e[c],r,t):-1===e.indexOf(r)&&a.push(cloneIfNecessary(r,t))}),a}function mergeObject(e,r,t){var a={};return isMergeableObject(e)&&Object.keys(e).forEach(function(r){a[r]=cloneIfNecessary(e[r],t)}),Object.keys(r).forEach(function(c){isMergeableObject(r[c])&&e[c]?a[c]=deepmerge(e[c],r[c],t):a[c]=cloneIfNecessary(r[c],t)}),a}function deepmerge(e,r,t){var a=Array.isArray(r),c=(t||{arrayMerge:defaultArrayMerge}).arrayMerge||defaultArrayMerge;return a?Array.isArray(e)?c(e,r,t):cloneIfNecessary(r,t):mergeObject(e,r,t)}deepmerge.all=function(e,r){if(!Array.isArray(e)||e.length<2)throw new Error("first argument should be an array with at least two elements");return e.reduce(function(e,t){return deepmerge(e,t,r)})};

/*

* ChemSpace.js is an easy-to-use versatile tool for the 2D visualization
* of compound sets within a web page.
* Source code, tutorial, documentation, and example
* data are freely available from ChemSpace.js website <a
* href="http://openscreen.cz/software/ChemSpace"
* target=blank>http://openscreen.cz/software/ChemSpace</a>. At the
* website, you can also find a Python script <a
* href="http://openscreen.cz/software/chemspace/chemspacepy"
* target=blank>chemspace.py</a> which process and prepares <a href="http://openscreen.cz/software/ChemSpace/input_format"
* target=blank>input data for ChemSpace.js</a>.
*
* @author <a href="mailto:ctibor.skuta@img.cas.cz">Ctibor Škuta</a>
* @author <a href="mailto:petr.bartunek@img.cas.cz">Petr Bartůněk</a>
* @author <a href="mailto:svozild@vscht.cz">Daniel Svozil</a>
* @version dev
* @category 1
* @license ChemSpace.js http://openscreen.cz/software/ChemSpace Copyright 2015, Ctibor Škuta, Petr Bartůněk, Daniel Svozil Licensed under the MIT license.
*
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
* OTHER DEALINGS IN THE SOFTWARE.
*
* @requires <a href='https://code.jquery.com/jquery-3.3.1.min.js'>jQuery Core 3.3.1</a>
* @dependency <script language="JavaScript" type="text/javascript" src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
*
* @requires <a href='http://konvajs.github.io'>KonvaJS 1.7.6</a>
* @dependency <script language="JavaScript" type="text/javascript" src="https://cdn.rawgit.com/konvajs/konva/1.7.6/konva.min.js"></script>
*
* @param {Object} options An object with the options for the ChemSpace.js component.
*
* @option {string} target
*   identifier of the DIV tag where the component should be displayed
*/

var ChemSpace;
var _date = new Date();

(function($){
  ChemSpace = function(settings){
      var self = this;
      self.user_settings = settings;
      self.target_element = $("#" + settings.target);
      var target_width = self.target_element.width();

      /**
      * Default color scales
      * @name ChemSpace#colors
      */

      self.colors = {
          RdYlBu: ['rgb(215,25,28)', 'rgb(255,255,178)', 'rgb(44,123,182)'],
          RdYlGr: ['rgb(215,25,28)', 'rgb(255,255,178)', 'rgb(26,150,65)'],
          BuWhRd: ['rgb(33,113,181)', 'rgb(255,255,255)', 'rgb(215,25,28)'],
          RdLrBu: ['rgb(215,25,28)', 'rgb(254,229,217)', 'rgb(44,123,182)'],
          RdLrGr: ['rgb(215,25,28)', 'rgb(254,229,217)', 'rgb(35,139,69)'],
      }

      self._color_dicts = {};
      
      Object.keys(self.colors).forEach(function(key){
        self.colors[key].forEach(function(rgb){
            self._color_dicts[rgb] = self._parse_color(rgb);
        })
      });

      if(self.user_settings.color !== undefined){
        
        if(typeof(self.user_settings.color.scale) === "string" && self.colors[self.user_settings.color.scale] !== undefined){
            var cs = self.user_settings.color.scale;
            self.user_settings.color.scale = {"min": self.colors[cs][0], "middle": self.colors[cs][1], "max": self.colors[cs][2]};
        }
        else if(typeof(self.user_settings.color.scale) === "object"){
            ["min", "middle", "max"].forEach(function(param){
                const color = self._parse_color(self.user_settings.color.scale[param]);
                if(color !== undefined){
                    const rgb = "rgb(" + color.r + "," + color.g + "," + color.b + ")";
                    self._color_dicts[rgb] = color;
                    
                }
            });
        }
      }

      /**
      * Default values for the settings
      * @name ChemSpace#settings
      */

      self.settings = {
          "target" : "YourOwnDivId", //the ID of the target element
          "height" : 800, // the height of visualization
          "width" : target_width, // the width of visualization (default is the width of target element)
          "font": {
            "fill": "#333333",
            "fontStyle": "bold",
            "fontFamily": "Mono, sans-serif",
            "fontSize": 14,
            "lineHeight": 1.2
          }, // font
          "categories": {
              "draw": true,
              "order": false,
              "off_opacity": 0.3,
              "legend": {
                "fontFamily": "Mono, sans-serif",
                "fill": "#333333",
                "fontWeight": "bold",
                "fontSize": 14,
                "width": 200
              },
          },
          "shapes": { // default color and radius of points, the width of stroke
            "order": ["circle", "triangle", "square", "rhombus", "triangle-down"],
            "default": {
              "radius": 4,
              "strokeWidth": 1,
              "fill": "#C2C2C2",
              "stroke": null,
              "shadowForStrokeEnabled": false,
              "preventDefault": false
            },
            "circle":{
              "transformsEnabled": "position",
            },
            "square":{
              "sides": 4,
              "rotation": 45,
              "transformsEnabled": "all",
            },
            "triangle":{
              "sides": 3,
              "transformsEnabled": "position",
            },
            "rhombus":{
              "sides": 4,
              "transformsEnabled": "position",
            },
            "triangle-down":{
              "sides": 3,
              "transformsEnabled": "all",
              "rotation": 180
            },
          },
          "path": { // default color and radius of points, the width of stroke
            "strokeWidth": 1,
            "stroke": 'gray',
          },
          "links": { // default color and radius of points, the width of stroke
            "draw": true,
            "attrs": {
                "strokeWidth": 1.5,
                "stroke": '#777777',
            }
          },
          "link_width": { // color scale settings
            "scale": {"min": 2, "middle": 4, "max": 6},
            "value_type": "percentile", // color by percentile/value
            "params": {"min": 5, "max": 95, "middle": 50}
          },
          "color": { // color scale settings
            "scale": {"min": self.colors.RdLrGr[0], "middle": self.colors.RdLrGr[1], "max": self.colors.RdLrGr[2]}, // color scale
            "value_type": "percentile", // color by percentile/value
            "index": 2, // feature index/category
            "params": {"min": 5, "max": 95, "middle": 50}
          },
          "point_size": { // color scale settings
            "scale": {"min": 2, "middle": 4, "max": 6},
            "value_type": "percentile", // color by percentile/value
            "index": 3, // feature index/category
            "params": {"min": 5, "max": 95, "middle": 50}
          },
          "coordinates": { // coordinates, x, y feature indexes
            "x": 0,
            "y": 1,
          },
          "highlight_color": "black",
          "navigation_toggle": { // show/hide navigation elements
            "header": true,
            "axis_labels": true,
            "axis": true,
            "diagonal": true,
            "export_button": true,
            "color": true,
            "point_size": true,
            "point_count": true,
            "link_width": true
          },          
          "compounds": { // compound settings
            draw: true, // whether to draw compounds
            limit: 50, // maximum limit of drawn compounds
            size: 200, // default size of molecules
            tooltip_compound_size: 100, // size of compound in the tooltip
            smilesDrawer: {
              bondThickness: 1,
              fontSizeLarge: 7,
              fontSizeSmall: 4,
            }
          },
          "resolution": 0,
          "colors": ["#d62728", "#2ca02c", "#1f77b4", "#ff7f0e", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"]
      };
      

      self.update_settings(settings);
      
      if(settings.categories === undefined || settings.categories.legend === undefined || settings.categories.legend.fontFamily === undefined){
        self.settings.categories.legend.fontFamily = self.settings.font.fontFamily;
      }

      if(self.settings.shapes.default.strokeWidth === 0){
        self.settings.shapes.default.stroke = null;
      }

      for(var i = 0, len=self.settings.shapes.order.length; i<len; i++){
        var key = self.settings.shapes.order[i];
        self.settings.shapes[key] = deepmerge(self.settings.shapes[key], self.settings.shapes.default);
      };

      self.use_cache = true;
      self.target_element.css({"position": "relative", "font-family": self.settings.font.fontFamily, "color": self.settings.font.fill});
      self._prepare_objects_ref();

      self._settings = {
        "point_size": {
            "menu_shape": $("<div></div>").css({           
                "background-color": "#999",
            }),
            "order": 3,
            "label": "Point radius",
            "property_label":  "Radius",
        },
        "link_width": {
            "menu_shape": $("<div></div>").css({           
                "width": "15px",
                "border-bottom-style": "solid",
                "border-bottom-color": "#999"
            }),
            "order": 2,
            "adjusted_attr": "border-color",
            "label": "Edge width",
            "property_label": "Width",
        },
        "color": {
            "menu_shape": $("<div></div>").css({           
                "width": "100%",
                "border": "solid #D2D2D2 1px",
                "height": 15,
                "min-width": 100
            }),
            "order": 1,
            "adjusted_attr": "border-color",
            "label": "Color scale",
            "property_label": "Color",
        },
      }

      /**
      * Default function definitions for the ChemSpace events
      * @name ChemSpace#events
      */
      self.events = {
          /**
          * @name ChemSpace#point_click
          * @event
          * @param {function} function() callback function for click on a single point event
          * @eventData {string} string point ID
          * @eventData {object} event event object

          * @example
          * instance.events.point_click = (
          *    function(point_id, evt) {
          *       alert(point_id);
          *    }
          * );
          *
          */
          "point_click": function(point_id, evt){
            return;
          },

          /**
          * @name ChemSpace#path_click
          * @event
          * @param {function} function() callback function for click on a path event
          * @eventData {object} object path object representation
          * @eventData {object} event event object

          * @example
          * instance.events.path_click = (
          *    function(path, evt) {
          *       alert(path);
          *    }
          * );
          *
          */
          "path_click": function(path, evt){
            return;
          },

          /**
          * @name ChemSpace#link_click
          * @event
          * @param {function} function() callback function for click on a link event
          * @eventData {object} object point ids linked by the link
          * @eventData {object} event event object

          * @example
          * instance.events.link_click = (
          *    function(point_ids, evt) {
          *       alert(point_ids);
          *    }
          * );
          *
          */
          "link_click": function(point_ids, evt){
            return;
          },

          /**
          * @name ChemSpace#points_selection
          * @event
          * @param {function} function() callback function for selection of points event
          * @eventData {array} array array of point IDs
          * @eventData {object} event event object

          * @example
          * instance.events.points_selection = (
          *    function(point_ids) {
          *       alert(point_ids);
          *    }
          * );
          *
          */
          "points_selection": function(point_ids){
            return;
          },

          /**
          * @name ChemSpace#category_legend_click
          * @event
          * @param {function} function() callback function for click on a category legend event
          * @eventData {object} object category object representation
          * @eventData {object} event event object

          * @example
          * instance.events.category_legend_click = (
          *    function(category_object, evt) {
          *       alert(category_object);
          *    }
          * );
          *
          */
          "category_legend_click": function(category_object, evt){
            return;
          },

          /**
          * @name ChemSpace#refresh
          * @event
          * @param {function} function() callback function for click on a refresh icon
          * @eventData {array} point ids
          * @eventData {object} event event object

          * @example
          * instance.events.refresh = (
          *    function(ids, evt) {
          *       alert(ids.length);
          *    }
          * );
          *
          */
          "refresh": function(category_object, evt){
            return;
          },

          /**
          * @name ChemSpace#point_tooltip
          * @event
          * @param {function} function() callback function for point tooltip customization
          * @eventData {object} event event object

          * @example
          * instance.events.point_tooltip = (
          *    function(evt) {
          *       alert(evt);
          *    }
          * );
          *
          */
          "point_tooltip": function(point_ids, color, evt){
            return self._get_point_tooltip(evt);
          },

      }

      self.paths_ref = {
            "zoom_icon": "M22.646,19.307c0.96-1.583,1.523-3.435,1.524-5.421C24.169,8.093,19.478,3.401,13.688,3.399C7.897,3.401,3.204,8.093,3.204,13.885c0,5.789,4.693,10.481,10.484,10.481c1.987,0,3.839-0.563,5.422-1.523l7.128,7.127l3.535-3.537L22.646,19.307zM13.688,20.369c-3.582-0.008-6.478-2.904-6.484-6.484c0.006-3.582,2.903-6.478,6.484-6.486c3.579,0.008,6.478,2.904,6.484,6.486C20.165,17.465,17.267,20.361,13.688,20.369zM15.687,9.051h-4v2.833H8.854v4.001h2.833v2.833h4v-2.834h2.832v-3.999h-2.833V9.051z",
            "unzoom_icon": "M22.646,19.307c0.96-1.583,1.523-3.435,1.524-5.421C24.169,8.093,19.478,3.401,13.688,3.399C7.897,3.401,3.204,8.093,3.204,13.885c0,5.789,4.693,10.481,10.484,10.481c1.987,0,3.839-0.563,5.422-1.523l7.128,7.127l3.535-3.537L22.646,19.307zM13.688,20.369c-3.582-0.008-6.478-2.904-6.484-6.484c0.006-3.582,2.903-6.478,6.484-6.486c3.579,0.008,6.478,2.904,6.484,6.486C20.165,17.465,17.267,20.361,13.688,20.369zM8.854,11.884v4.001l9.665-0.001v-3.999L8.854,11.884z",
            "options_icon": "M4.082,4.083v2.999h22.835V4.083H4.082zM4.082,20.306h22.835v-2.999H4.082V20.306zM4.082,13.694h22.835v-2.999H4.082V13.694zM4.082,26.917h22.835v-2.999H4.082V26.917z",
            "refresh_icon": "M24.083,15.5c-0.009,4.739-3.844,8.574-8.583,8.583c-4.741-0.009-8.577-3.844-8.585-8.583c0.008-4.741,3.844-8.577,8.585-8.585c1.913,0,3.665,0.629,5.09,1.686l-1.782,1.783l8.429,2.256l-2.26-8.427l-1.89,1.89c-2.072-1.677-4.717-2.688-7.587-2.688C8.826,3.418,3.418,8.826,3.416,15.5C3.418,22.175,8.826,27.583,15.5,27.583S27.583,22.175,27.583,15.5H24.083z",
            "export_icon": "M24.25,10.25H20.5v-1.5h-9.375v1.5h-3.75c-1.104,0-2,0.896-2,2v10.375c0,1.104,0.896,2,2,2H24.25c1.104,0,2-0.896,2-2V12.25C26.25,11.146,25.354,10.25,24.25,10.25zM15.812,23.499c-3.342,0-6.06-2.719-6.06-6.061c0-3.342,2.718-6.062,6.06-6.062s6.062,2.72,6.062,6.062C21.874,20.78,19.153,23.499,15.812,23.499zM15.812,13.375c-2.244,0-4.062,1.819-4.062,4.062c0,2.244,1.819,4.062,4.062,4.062c2.244,0,4.062-1.818,4.062-4.062C19.875,15.194,18.057,13.375,15.812,13.375z",
      };

      self.shape_events = {
        "mouseover":{
          "point": function(evt){self._point_mouseover(evt);},
          "path": function(evt){self._path_mouseover(evt);},
          "link": function(evt){self._link_mouseover(evt);},
          "Rect": function(evt){return},
          "category_legend": function(evt){self._category_legend_mouseover(evt);},
          "icon": function(evt){self._icon_mouseover(evt);},
        },
        "mouseout":{
          "point": function(evt){self._point_mouseout(evt);},
          "path": function(evt){self._point_mouseout(evt);},
          "link": function(evt){self._link_mouseout(evt);},
          "category_legend": function(evt){self._category_legend_mouseout(evt);},
          "icon": function(evt){self._icon_mouseout(evt);},
          "Rect": function(evt){return},
        },
        "click":{
          "point": function(evt){self._point_click(evt);},
          "path": function(evt){self._path_click(evt);},
          "link": function(evt){self._link_click(evt);},
          "category_legend": function(evt){self._category_legend_click(evt);},
          "icon": function(evt){self._icon_click(evt);},
          "Rect": function(evt){return},
        }
      };

      /*self.shape_events.touchstart = self.shape_events.mouseover;
      self.shape_events.touchend = self.shape_events.mouseout;*/

      self.html_ref = {
        "tooltip": $("<div class='point_tooltip'>\
            <div class='compound_img'><canvas></canvas></div>\
            <div class='point_id'></div>\
            <div class='dimensions'></div>\
          </div>")
          .css({
            "border": "solid #D2D2D2 3px",
            "padding": 10,
            "font-size": "80%",
            "background-color": "white",
            "word-break": "break-all",
            "max-width": self.settings.compounds.tooltip_compound_size+20
          }),

        "redraw_button": $('<div class="redraw_button">Redraw</div>')
          .css({
            "padding": "0.5rem",
            "color": "white",
            "border-radius": "2px",
            "margin": "0.5rem 0rem",
            "background-color": "#333",
            "font-weight": "bold",
            "text-align": "center",
          }),

          "menu": $("<div></div>")
            .css({"position": "absolute",
              "border": "solid #D2D2D2 2px",
              "padding": "4px 10px",
              "font-size": "small",
              "background-color": "white",
              "z-index": 30,
              "color": "#333333"
            }),

          "menu_subtitle": $("<div></div>").css({"font-weight": "bold", "padding": "5px 0px"}),
      }
      
      self.tooltip_templates = {
        "default": function(point_ids, tooltip_wrapper){
            for(var i1 = 0, len1=point_ids.length; i1<len1; i1++){
                var objects = self.id2points[point_ids[i1]];
                var coords = self.point2coord[objects[0]];
                var object_ids = self.point_index[coords[0]][coords[1]];
                var object_id = objects[0];
                var label = self.data.points[object_id][self.keys.label];
                label = $("<div><div>" + ((label === undefined)?object_id:label) + "</div></div>").css({"display": "flex", "align-items": "center"});
                var point_data = self._get_points_data(object_ids);     
                
                if(object_ids.length > 1){
                  label.append($("<div>+ " + (object_ids.length-1) + "</div>")
                    .css({"border-radius": 5, "padding": "3px 8px", "background-color": "#65aa30", "color": "white", "margin-left": 8, "white-space": "nowrap"})
                  );
                }
                
                var color = self.stage.find("#"+point_ids[i1])[0].fill();
                var tooltip = self.html_ref.tooltip.clone();
                tooltip.find("canvas").attr("id", self.settings.target + "@tooltip_" + i1);

                tooltip.find(".point_id").html(label)
                  .css({
                    "font-weight": "bold",
                    "border-bottom": "solid #d2d2d2 1px",
                    "padding-bottom": 5,
                });

                tooltip.css({
                  "border-color":color,            
                  "margin-left": ((i1 !== 0)?5:0)
                });

                var dimensions = $("<div></div>");

                for(var i2 = 0, len2=self.data.feature_names.length; i2<len2; i2++){
                  row = $("<div></div>")
                    .css({"margin-top": 5});

                  row.append($("<span>" + self.data.feature_names[i2] + ": </span>: ").css({"color": "#333333", "margin-right": 5}));
                  row.append("<span>" + ((point_data[i2] !== null)?point_data[i2]:"-") + "</span>");

                  dimensions.append(row);
                }
                tooltip.find(".dimensions").html(dimensions);
                tooltip_wrapper.append(tooltip);
            }
            return tooltip_wrapper;
        }
      };

      self.target_element.append(self.html_ref.tooltip);

      if(self.settings.compounds.draw){

        self.img_cache = $("#chemspacejs-img_cache");
        if(self.img_cache.length == 0){
            self.img_cache = $("<canvas id='chemspacejs-img_cache'></canvas>").css({"position": "fixed", "top": -1000});
        }
        $("body").append(self.img_cache);

        var sd_space = $.extend({}, self.settings.compounds.smilesDrawer);
        sd_space.width = self.settings.compounds.size;
        sd_space.height = self.settings.compounds.size;
        self.smilesDrawer = new SmilesDrawer.Drawer(sd_space);

        var sd_tooltip = $.extend({}, self.settings.compounds.smilesDrawer);
        sd_tooltip.width = self.settings.compounds.tooltip_compound_size;
        sd_tooltip.height = self.settings.compounds.tooltip_compound_size;
        self.tooltipSmilesDrawer = new SmilesDrawer.Drawer(sd_tooltip);
      
      //   self.compound_img_stage = new Konva.Stage({
      //     container: self.target_element.find(".compound_img")[0],
      //     width: self.settings.compounds.tooltip_compound_size+10,
      //     height: self.settings.compounds.tooltip_compound_size+10
      //   });
      //   self.compound_img_layer = new Konva.Layer({"listening": false});
      //   self.compound_img_stage.add(self.compound_img_layer);
      }

      self._prop2settings = {
        "color": $.extend(true, {}, self.settings.color),
        "point_size": $.extend(true, {}, self.settings.point_size),
        "link_width": $.extend(true, {}, self.settings.link_width),
        "coordinates": $.extend(true, {}, self.settings.coordinates),
      };

      self._prop2settings.color.scale = "init";
      self._prop2settings.point_size.scale = "init";
      self._prop2settings.link_width.scale = "init";
      self._prop2settings.coordinates.x = -1;

      self._prop2fnc = {
        "color": {
          "category": function(category){
            return category.color;
          },
          "feature": function(value){
            return self._get_color_for_value(
                value,
                self._prop2settings.color.values.min,
                self._prop2settings.color.values.max,
                self._prop2settings.color.values.middle,
                self.settings.color.scale);
          }
        },
        "point_size":{
          "category": function(category){
            return category.radius;
          },
          "feature": function(value){
            return self._get_value_on_scale(value, "point_size");
          }
        },
        "link_width":{
          "feature": function(value){
            return self._get_value_on_scale(value, "link_width");
          }
        }
      };
  }

  ChemSpace.prototype._prepare_objects_ref = function(){
    /**
      * Default konvajs objects references
      * @name ChemSpace#objects_ref
      */
      var self = this;
      self.shapes_ref = {};
      
      for(var i = 0, len=self.settings.shapes.order.length; i<len; i++){
        var key = self.settings.shapes.order[i];
        var shape_settings = $.extend(self.settings.shapes.default, {class: "point"});
        shape_settings = deepmerge(shape_settings, self.settings.shapes[key]);

        if(key === 'circle'){
          self.shapes_ref[key] = new Konva.Circle(shape_settings);
        }
        else{
          self.shapes_ref[key] = new Konva.RegularPolygon(shape_settings);
        }
      };

      self.keys = {
        "default": {
          "object_ids": "object_ids",
          "label": "label",
          "smiles": "smiles",
          "features": "features",
        },
        "compressed":{
          "object_ids": "o",
          "label": "l",
          "smiles": "s",
          "features": "f",
        }
      }

      self.objects_ref = {
          "path": new Konva.Line({
            stroke: 'gray',
            strokeWidth: self.settings.path.strokeWidth,
            lineCap: 'round',
            lineJoin: 'round',
            opacity: 0.7,
            class: "path"
          }),

          "link": new Konva.Line({
            stroke: self.settings.links.attrs.stroke,
            strokeWidth: self.settings.links.attrs.strokeWidth,
            lineCap: 'round',
            lineJoin: 'round',
            opacity: 0.7,
          }),

          "link_hit_area": new Konva.Line({
            stroke: "#FFFFFF",
            opacity: 0,
            class: "link"
          }),

          "selection_rect": new Konva.Rect({
            width: 0,
            height: 0,
            fill: "#D2D2D2",
            opacity: 0.3,
            listening: false,
          }),

          "icon": new Konva.Path({
            fill: "grey",
            class: "icon",
            listening: false
          }),

          "icon_overlay": new Konva.Rect({
            width: 32,
            height: 32,
            opacity: 0,
            class: "icon",
          }),

          "legend_text": new Konva.Text(self.settings.font),

          "legend_circle": new Konva.Circle({
            radius: 6,
            strokeWidth: 3
          }),

          // "hover_overlay_rect": new Konva.Rect({
          //   listening: false,
          //   opacity: 0.8,
          //   fill: "white",
          //   width: self.settings.width,
          //   x: self.left_margin,
          // }),

          "tooltip_label": new Konva.Label({
            opacity: 1,
            listening: false,
          }),

          "tooltip_tag": new Konva.Tag({
            fill: "gray",
            pointerWidth: 10,
            pointerHeight: 10,
            lineJoin: 'round',
            listening: false,
          }),

          "tooltip_text": new Konva.Text({
            fontFamily: self.settings.font.fontFamily,
            fontSize: self.settings.font.fontSize,
            padding: 8,
            fill: 'white',
            fontStyle: "bold",
            listening: false,
            align: "center",
            lineHeight: 1.2,
          }),

          "arrow": new Konva.Arrow({
            pointerLength: 10,
            pointerWidth : 7,
            fill: '#B2B2B2',
            stroke: '#B2B2B2',
            listening: false
          }),

          "diagonal": new Konva.Line({
            stroke: "#B2B2B2",
            strokeWidth: 2,
            lineCap: "round",
            dash: [10,5],
            listening: false
          }),

          "border_line": new Konva.Line({
            stroke: "#B2B2B2",
            strokeWidth: 1,
            lineCap: "round",
            dash: [4,2],
            listening: false
          }),

          "rect_gradient": new Konva.Rect({
            x: 0,
            y: 80,
            height: 20,
            fillLinearGradientStartPoint: {x: 0, y: 80},
            stroke: "#333333",
            strokeWidth: 0.5
          }),

          "percentile_line": new Konva.Line({
            stroke: "#393939",
            strokeWidth: 0.5,
            lineCap: "round"
          }),

          "navigation_text": new Konva.Text(self.settings.font),
      };

      self.objects_ref.shapes = self.shapes_ref;

      self.shape_fnc = function(shape, point_id, point_ids, x, y, radius, color){
          return self.objects_ref.shapes[shape].clone({
              x: x,
              y: y,
              point_ids: point_ids,
              fill: color,
              radius: radius,
              id: point_id
          });
      };
  }

  /**
    * Update ChemSpace settings.
    *
    * @param {Object} [settings] Settings object in the same format as for the initialization.
    *
    */
  ChemSpace.prototype.update_settings = function(settings){
    var self = this;
    // if((settings.width !== undefined && self.settings.width !== settings.width) || (settings.height !== undefined && self.settings.height !== settings.height)){
        self.use_cache = false;
    // }
    self.settings = deepmerge(self.settings, settings);
  }

  /**
    * Read data from JSON variable.
    *
    * @param {object} [variable] Chemical space in proper JSON format.
    */
  ChemSpace.prototype.read_data = function(json){
      var self = this;
      self.data = json;
      console.log(self.data)

      self.point_ids = Object.keys(self.data.points);
      self.points_len = self.point_ids.length;

      if(self.points_len > 0){
        if(self.data.points[0]["o"] !== undefined){
          self.keys = self.keys.compressed;
        }
        else{
          self.keys = self.keys.default;
        }

      }

      self.object_id2points = {};

      for(var i = 0, len=self.points_len; i<len; i++){
        var point_id = self.point_ids[i];
        for(var j = 0, len_j=self.data.points[point_id][self.keys.object_ids].length; j<len_j; j++){
          var object_id=self.data.points[point_id][self.keys.object_ids][j];
          if(self.object_id2points[object_id] === undefined){
            self.object_id2points[object_id] = [];
          }
          self.object_id2points[object_id].push(point_id);
        }
      }
      self._process_categories(self.data.categories);
      self.point_id2categories = self._get_id2categories();

      self._get_links_feature();

      if(self.settings.compounds.draw && self.fetch_smiles === undefined){
          self.settings.compounds.draw = false;
      }

      var point = self.data.points[self.point_ids[0]];
      var dims = point[self.keys.features].length;

      if(self.data.feature_names === undefined){
        self.data.feature_names = [];
        for(var i = 0; i<dims; i++){
          self.data.feature_names.push("Feature " + (i+1));
        }
      }

      // self.objects_ref.hover_overlay_rect.setAttrs({y: self.header_height, height: self.settings.height - self.header_height});
  }

  /**
    * Read data from JSON file.
    *
    * @param {string} [filename] Path to the JSON data file.
    *
    */
  ChemSpace.prototype.read_data_from_file = function(json){
    var self = this;
      $.ajax({
          type: 'GET',
          url: json,
          dataType: 'json',
          success: function(json_file){
            self.read_data(json_file);
          },
          async: false
      });
  }

  /*ChemSpace.prototype.fetch_smiles = function(){
    return false;
  }*/  

  ChemSpace.prototype._add_prefix_to_data = function(data){
    var self = this;
    var id, prefixed_data = {};

    for(var i = 0, keys = Object.keys(data), len = keys.length; i < len; i++){
        id = [self.settings.target, keys[i]].join("#");
        prefixed_data[id] = data[keys[i]];
    }
    return prefixed_data;
  }

  ChemSpace.prototype._draw_stage_layer = function(){
    var self = this;
      self.stage_layer = new Konva.Layer();
      self.stage_rect = new Konva.Rect({
                                  x: 0,
                                  y: 0,
                                  width: self.settings.width,
                                  height: self.settings.height,
                                  opacity: 0,
                                  preventDefault: false
                              });
      self.stage_layer.add(self.stage_rect);
      self.stage_rect.moveToBottom();
      self.stage.add(self.stage_layer);

      self.stage_rect.on("click", function(evt){
        if(self.selection_rect){
          self.selection_rect.destroy();
          self.zoom_icon.destroy();
          self.selection_layer.draw();
        }
      });
  }

  ChemSpace.prototype._destroy_children = function(layers){
    var self = this;
    for(var i = 0, len=layers.length; i<len; i++){
      layers[i].destroyChildren();
    }
  }

   /**
    * Add feature.
    *
    * @param {Object} [category] Category defined by label and points in format {label: string, "color": string, "points": [point_id, point_id ...]}.
    *
    */
  ChemSpace.prototype.add_feature = function(feature){
    var self = this;
    if(self.data.feature_names.indexOf(feature.name) > -1){
      alert("Feature '" + feature.name + "' already exists.");
    }
    else{
        for(var i = 0, len = self.points_len; i < len; i++){
            var point_id = self.point_ids[i];
            var object_ids = self.data.points[point_id][self.keys.object_ids];
            var values = [];
            var total = 0;
            
            for(var object_i = 0, object_len = object_ids.length; object_i < object_len; object_i++){
              if(feature.point2value[object_ids[object_i]] !== undefined){
                  values.push(feature.point2value[object_ids[object_i]]);
                  total += feature.point2value[object_ids[object_i]];
              }
            }

            if(values.length == 1){
                self.data.points[point_id][self.keys.features].push(values[0]);
            }
            else if(values.length > 1){
                self.data.points[point_id][self.keys.features].push(self._hack_round(total/values.length));
            }
            else{
                self.data.points[point_id][self.keys.features].push(null);
            }
        }
        self.data.feature_names.push(feature.name);
        var selects = self.target_element.find(".navigation select");
        var feature_index = self.data.feature_names.length-1;

        for(var i = 0, len=selects.length; i<len; i++){
            $(selects[i]).append("<option value='" + feature_index + "''>" + feature.name + "</option>");
        }
    }
  }

   /**
    * Remove feature.
    *
    * @param {String} [name] Name of the feature to be removed.
    *
    */
  ChemSpace.prototype.remove_feature = function(name){
    var self = this;
    var f_index = self.data.feature_names.indexOf(name);
    if(f_index > -1){
      self.data.feature_names.splice(f_index, 1);

      for(var i = 0, len = self.points_len; i < len; i++){
        var point_id = self.point_ids[i];
        self.data.points[point_id][self.keys.features].splice(f_index, 1);
      }

      self.target_element.find(".navigation select option[value='" + f_index + "']").remove();
    }

  }

   /**
    * Add category.
    *
    * @param {Object} [category] Category defined by label and points in format {label: string, "color": string, "points": [point_id, point_id ...], "objects": [object_id, object_id ...]}.
    *
    */
  ChemSpace.prototype.add_category = function(category){
    var self = this;
    if(self.data.categories === undefined){
      self.data.categories = [];
    }
    self.categories = [];
    self.category2layer = {};
    self.data.categories.push(category);
    self._process_categories(self.data.categories);
    self.point_id2categories = self._get_id2categories();
    // self._sort_layers();
  }

  /**
    * Update category.
    *
    * @param {Object} [category] Category defined by label and points in format {label: string, "color": string, "points": [point_id, point_id ...]}.
    *
    */
  ChemSpace.prototype.update_category = function(category, action){
    var self = this;
    if(action === undefined){
      var action = "update";
    }
    var update = false;

    self._prop2settings.color.scale = "init";
    self._prop2settings.point_size.scale = "init";
    self._prop2settings.coordinates.x = -1;

    if(self.data.categories !== undefined){
      var len = self.data.categories.length;

      for(var i = 0; i<len; i++){
        var c = self.data.categories[i];

        if(c.label === category.label){
          if(action === "update"){
            self.data.categories[i] = category;
          }
          else if(action === "remove"){
            self.data.categories.splice(i, 1);
          }
          update = true;
          break;
        }
      }

      if(update){

        self.categories = [];
        self.category2layer = {};
        self._process_categories(self.data.categories);
        self.point_id2categories = self._get_id2categories();
        self._sort_layers();
      }
      else{
        self.add_category(category);
      }
    }
    else{
      self.add_category(category);
    }
  }

  /**
    * Remove category.
    *
    * @param {String} [label] Remove category by label.
    *
    */
  ChemSpace.prototype.remove_category = function(label){
    var self = this;
    self.update_category({label: label}, "remove");
  }

   /**
    * Add path.
    *
    * @param {Object} [path] Path defined by points in format {label: string, "color": string, "points": [point_id, point_id ...]}.
    *
    */
  ChemSpace.prototype.add_path = function(path){
    var self = this;
    if(self.data.paths === undefined){
      self.data.paths = [];
    }
    self.data.paths.push(path);
    if(self.path_layer != undefined){
      self.path_layer.destroyChildren();
    }
    else{
      self.path_layer = new Konva.Layer();
      self.stage.add(self.path_layer);
    }
    self._draw_paths();
    self._sort_layers();
  }

  /**
    * Unhighlight points.
    * Unhighlights all highlighted points.
    *
    * @example
    * instance.unhighlight_points();
    */

  ChemSpace.prototype.unhighlight_points = function(){
    var self = this;
    self.highlight_layer.destroyChildren();
    self.highlight_layer.draw();
  }

  /**
    * Highlight points.
    * When the empty array is passed it unhighlights all highlighted points.
    *
    * @param {object} [point_ids] The array of point (object) IDs.
    *
    * @example
    * instance.highlight_rows(["id1", "id2"]);
    */

  ChemSpace.prototype.highlight_points = function(point_ids, color){
    var self = this;
    var len = point_ids.length;
    if(color !== undefined){
      self.settings.highlight_color = color;
    }

    if(self.highlight_layer === undefined){
      self.highlight_layer = new Konva.Layer();
      self.stage.add(self.highlight_layer);
    }
    else{
      self.highlight_layer.destroyChildren();
    }

    self.highlighted_points = point_ids;

    if(len > 0){
      for(var i = 0; i<len; i++){
        self.highlight_layer.add(self.stage.find("#"+self.point2id[self.highlighted_points[i]])[0].clone({"listening": false, "fill": self.settings.highlight_color}));
      }
      self.highlight_layer.draw();
    }
    else{
      self.unhighlight_points();
    }

  }

  /**
    * Redraw Chemspace.
    */
  ChemSpace.prototype.redraw = function(){
    var self = this;
    self.stage.destroy();
    self.stage = undefined;
    self._prepare_objects_ref();
    self.draw();
    self.use_cache = true;
  }

  ChemSpace.prototype.redraw_points = function(point_ids){
    var self = this;

    if(self.path_layer !== undefined){
      self.path_layer.destroyChildren();
    }
    self._draw_points(point_ids);

    if(self.data.paths !== undefined){
      self._draw_paths();
    }
    self._sort_layers();
  }

  ChemSpace.prototype._prepare_canvas = function(){
    var self = this;
    self.current_selection = [];

    self.stage = new Konva.Stage({
        container: self.settings.target,
    });

    self.stage.setWidth(self.settings.width);
    self.stage.setHeight(self.settings.height);
    self.target_element.height(self.settings.height);
    self.target_element.width(self.settings.width);

    self._draw_stage_layer();
  }

  ChemSpace.prototype.draw = function(){
    var self = this;
    if(self.stage === undefined){
      self._prepare_canvas();
    }
    self.header_height = (self.settings.navigation_toggle.header)?80:10;
    self.footer_height = (self.settings.navigation_toggle.axis_labels)?self.settings.font.fontSize + 10:10;
    self.left_margin = (self.settings.navigation_toggle.axis_labels)?self.settings.font.fontSize + 10:10;
    self.right_margin = 20;
    
    if(self.categories.length > 1 && self.settings.categories.draw){
      self.right_margin = self.right_margin + self.settings.categories.legend.width;
    }

    self.highlighted_points = [];
    self.fixed_compounds = [];
    self.selection_overlay_rect = new Konva.Rect({
      fill: "white",
      opacity: 0,
      x:self.left_margin,
      y: self.header_height,
      width: self.settings.width-self.right_margin-self.left_margin,
      height: self.settings.height-self.header_height-self.footer_height,
      preventDefault: false
    });
    self.current_selection = {
      "x": [self.left_margin, self.settings.width-self.right_margin-self.left_margin],
      "y": [self.header_height, self.settings.height-self.header_height-self.footer_height]
    }

    self.main_layer = new Konva.Layer();
    self.main_layer.add(self.selection_overlay_rect);
    self.link_layer = new Konva.Layer();
    self.compounds_layer = new Konva.FastLayer({"listening": false});
    self._draw_points();
    self._draw_menu();
    self._draw_legend();
    self.point2img = {};

    if(self.settings.compounds.draw){
      self._draw_compounds();
    }

    self.hover_layer = new Konva.Layer({"listening": false});
    self.stage.add(self.link_layer, self.main_layer, self.hover_layer, self.compounds_layer);

    if(self.data.paths !== undefined){
      self.path_layer = new Konva.Layer();
      self._draw_paths();
      self.stage.add(self.path_layer);
    }

    self._bind_hover_events_on_stage();

    self.selection_rect = false;
    self.down = false;

    self.selection_layer = new Konva.Layer();
    self.stage.add(self.selection_layer);
    self.selection_layer.moveToBottom();
    self.selection_background_rect = self.selection_overlay_rect.clone({opacity:0, listening: false});
    self.selection_layer.add(self.selection_background_rect);

    self.stage_layer.moveToBottom();
    self.stage.draw();
    self._sort_layers();

    self.target_element.find("input").css({"font-family": self.settings.font.fontFamily});
    self.target_element.css({"font-size": self.settings.font.fontSize});

    self.main_layer.on("mousedown", function(evt) {
      self.selection_background_rect.listening(true);
      self.selection_layer.draw();

      if(self.selection_rect){
        self.selection_rect.destroy();
        self.selection_rect = false;
        self.zoom_icon.destroy();
      }

      self.down = true;
      self.selection_layer.moveToTop();
      self.selection_rect = self.objects_ref.selection_rect.clone({x: evt.evt.offsetX, y: evt.evt.offsetY});
      self.selection_layer.add(self.selection_rect);

      self.stage.off("mouseover");
      self.stage.off("mouseout");
      self.stage.off("click");

      $(document).one("mouseup", function(){
          self.selection_layer.fire("mouseup");
      });
    });

    self.selection_layer.on("mousemove", function(evt) {
        if (!self.down) return;
        var evt = evt.evt;
        var rect = self.selection_rect.attrs;
        var width = evt.layerX - rect.x;
        var height = evt.layerY - rect.y;

        if(width*height !== 0){
          self.stage_layer.off("click");
        }

        self.selection_rect.width(width);
        self.selection_rect.height(height);
        self.selection_layer.draw();
    });

    self.selection_layer.on("mouseup", function(evt) {
        if(!self.down){
          return;
        }

        self.selection_background_rect.listening(false);
        self.selection_layer.draw();
        self.down = false;

        var rect = self.selection_rect.attrs;
        if(rect === undefined){
          return;
        }

        self.selection_range = {x:[rect.x, rect.x+rect.width], y:[rect.y, rect.y+rect.height]};
        self.selection_range.x.sort(self._sort_number_ascending);
        self.selection_range.y.sort(self._sort_number_ascending);
        self.current_selection = self.get_points_by_range(self.selection_range);
        self._bind_hover_events_on_stage();

        if(rect.width === 0||self.current_selection.length === 0){
          self.selection_rect.destroy();
          self.selection_rect = false;
          self.selection_layer.draw();
          return;
        }

        var rect_width = self.selection_rect.width();
        rect_width = (rect_width > 0)?rect_width:0;
        var rect_height = self.selection_rect.height();
        rect_height = (rect_height < 0)?rect_height:0;
        var x = self.selection_rect.getAttr("x") + rect_width + 5;
        var y = self.selection_rect.getAttr("y") + rect_height + 5;

        self.zoom_icon = self._get_icon_group("zoom_icon", "Zoom in", {x: x, y: y});
        self.zoom_icon.children[1].setAttrs({"fill": "white", "opacity": 0.7});
        self.zoom_icon.children[1].moveToBottom();
        self.selection_layer.add(self.zoom_icon);
        self.selection_layer.draw();

        self.events.points_selection(self._translate_points_to_objects(self.current_selection));

    });
  }

  ChemSpace.prototype._bind_hover_events_on_stage = function(){
    var self = this;

    self.stage.on("click", function(evt){
      if(evt.target.attrs.class !== undefined && self.shape_events[evt.type][evt.target.attrs.class] !== undefined){
        self.shape_events[evt.type][evt.target.attrs.class](evt);
      }
    });

    self.stage.on("mouseover", function(evt){
      if(evt.target.attrs.class !== undefined && self.shape_events[evt.type][evt.target.attrs.class] !== undefined){
        self.shape_events[evt.type][evt.target.attrs.class](evt);
      }
    });

    self.stage.on("mouseout", function(evt){
      if(evt.target.attrs.class !== undefined && self.shape_events[evt.type][evt.target.attrs.class] !== undefined){
        self.shape_events[evt.type][evt.target.attrs.class](evt);
      }
    });
  }

  ChemSpace.prototype._draw_points = function(point_ids){
    console.time("DRAW POINTS");
    var self = this;
    var keys, key, category, point, color, radius, point_ids,
      point_id, category_id, color_value, radius_value, categories,
      x_keys, x_key, y_keys, y_key, obj, category2data, point_count;

    var prop2cache = {
      coordinates: self._calculate_coordinates(point_ids),
      color: self._get_prop_settings("color"),
      // link_width: self._get_prop_settings("link_width"),
      point_size: self._get_prop_settings("point_size")
    };

    if(!prop2cache.coordinates || !prop2cache.point_size || !prop2cache.color){
      self._destroy_children(Object.values(self.category2layer));
      self._id2object = {};
      self.category2count = {};
      self.point2id = {};
      self.id2points = {};

      for(var i = 0, len=self.categories.length; i<len; i++){
        self.category2count[i] = 0;
      }

      x_keys = self._shuffle(Object.keys(self.point_index));
      var gen_id = 0;

      for(var i = 0, x_len=x_keys.length; i < x_len; i++){
        x_key = x_keys[i];
        y_keys = Object.keys(self.point_index[x_key]);

        for(var j = 0, y_len=y_keys.length; j < y_len; j++){
          y_key = y_keys[j];
          point_ids = self.point_index[x_key][y_key];
          point_count = point_ids.length;
          category2data = {};
          categories = [];

          for(var p = 0, p_len=point_count; p<p_len; p++){
            point_id = point_ids[p];            

            for(var c = 0, c_len=self.point_id2categories[point_id].length; c<c_len; c++){
              category_id = self.point_id2categories[point_id][c];

              if(category2data[category_id] === undefined){
                categories.push(category_id);
                category2data[category_id] = {
                  "color": ((self.settings.color.index != "category")?[]:category_id),
                  "radius": ((self.settings.point_size.index != "category")?[]:category_id),
                  "points": [],
                  "id": [self.settings.target, gen_id].join("@"),
                }
                gen_id++;
              }

              category2data[category_id].points.push(point_id);
              self.point2id[point_id] = category2data[category_id].id;

              if(self.settings.color.index != "category"){
                category2data[category_id].color.push(self.data.points[point_id][self.keys.features][self.settings.color.index]);
              }

              if(self.settings.point_size.index != "category"){
                category2data[category_id].radius.push(self.data.points[point_id][self.keys.features][self.settings.point_size.index]);
              }
            }
          }

          for(var c = 0, c_len=categories.length; c<c_len; c++){
            var cid = categories[c];
            if(self.settings.color.index != "category"){
              color = self._prop2fnc.color.feature(self._get_average(category2data[cid].color));
            }
            else{
              color = self._prop2settings.color.v2v[cid];
            }

            if(self.settings.point_size.index != "category"){
              radius = self._prop2fnc.point_size.feature(self._get_average(category2data[cid].radius));
            }
            else{
              radius = self._prop2settings.point_size.v2v[cid];
            }
            var shape = (self.categories[cid].shape !== undefined)?self.categories[cid].shape:"circle";

            self.id2points[category2data[cid].id] = category2data[cid].points;
            point = self.shape_fnc(shape, category2data[cid].id, category2data[cid].points, x_key, y_key, radius, color);
            self._id2object[category2data[cid].id] = point;
            self.category2count[cid] = self.category2count[cid] + category2data[cid].points.length;
            self.category2layer[cid].add(point);
          }


            // for(var c = 0, c_len=self.point_id2categories[point_id].length; c<c_len; c++){
            //   category_id = self.point_id2categories[point_id][c];
            //   category = self.categories[category_id];
            //   color_value = (self.settings.color.index != "category")?self.data.points[point_id].features[self.settings.color.index]:category_id;
            //   color = self._prop2settings.color.v2v[color_value];
            //   radius_value = (self.settings.point_size.index != "category")?self.data.points[point_id].features[self.settings.point_size.index]:category_id;
            //   radius = self._prop2settings.point_size.v2v[radius_value];

            //   point = self.shape2fnc.circle(point_id, x_key, y_key, radius, color);
            //   self._id2object[[point_id, category_id].join("_")] = point;
            //   self.category2count[category_id]++;
            //   self.category2layer[category_id].add(point);
            //   counter++;
            // }
          // }
        }
      }

      for(var i = 0, len=self.categories.length; i<len; i++){
        self.stage.add(self.category2layer[i]);
      };

      self._destroy_children([self.link_layer]);
      self._draw_links();
    }
    // else{
    //   if(!prop2cache.color){
    //     keys = Object.keys(self._id2object);
    //     for(var i = 0, len=keys.length; i<len; i++){
    //       key = keys[i];
    //       obj = self._id2object[key];
    //       category_id = parseInt(key.split("_")[1]);
    //       color_value = (self.settings.color.index != "category")?self.data.points[obj.attrs.id].features[self.settings.color.index]:category_id;
    //       obj.fill(self._prop2settings.color.v2v[color_value]);
    //       obj.draw();
    //     }
    //   }
    // }

    if(!prop2cache.coordinates){
      /*self._destroy_children([self.link_layer]);
      self._draw_links();*/
    }

    self.highlight_points(self.highlighted_points);

    console.timeEnd("DRAW POINTS");
  }

  ChemSpace.prototype._draw_links = function(point_ids){
    var self = this;

    if(!self.settings.links.draw){
      return;
    }

    var link, links, coords, p1, p2, hit_area, link_group, stroke_width;
    var weighted = false;
    if(point_ids === undefined){
      var point_ids = Object.keys(self.point2coord);
    }
    
    console.time("LINKS");
    for(var i = 0, len = point_ids.length; i < len; i++){
      p1 = point_ids[i];
      links = self.data.points[p1].links;

      if(links !== undefined && links.length > 0){
        if(typeof(links[0]) === "object"){
            self._get_prop_settings("link_width");
            weighted = true;
            break;
        }
        else{
            weighted = false;
            break;
        }
      }
    }

    for(var i = 0, len = point_ids.length; i < len; i++){
      p1 = point_ids[i];
      links = self.data.points[p1].links;

      if(links !== undefined && links.length > 0){
        if(weighted){
            for(var l = 0, len_links = links.length; l<len_links; l++){
                coords = self.point2coord[p1];
                p2 = links[l][0];
                if(self.point2coord[p2] !== undefined){
                    coords = coords.concat(self.point2coord[p2]);
                    link_group = new Konva.Group();
                    
                    stroke_width = self._prop2fnc.link_width.feature(links[l][1]);
                    link = self.objects_ref.link.clone({
                        "strokeWidth": stroke_width,
                        "points": coords,
                        "listening": false
                    });

                    hit_area = self.objects_ref.link_hit_area.clone({
                        "points": coords,
                        "strokeWidth": (stroke_width > 6)?stroke_width:6,
                        "point_ids": [self.point2id[p1], self.point2id[p2]],
                        "link": link,
                        "weight": links[l][1]
                    });

                    link_group.add(hit_area, link)
                    self.link_layer.add(link_group);
                  }
            }
        }
        else{
          for(var l = 0, len_links = links.length; l<len_links; l++){
              p2 = links[l];
              
              if(self.point2coord[p2] !== undefined){
                coords = self.point2coord[p1];
                coords = coords.concat(self.point2coord[p2]);
                link_group = new Konva.Group();
                link = self.objects_ref.link.clone({
                    "points": coords,
                    "listening": false
                });

                hit_area = self.objects_ref.link_hit_area.clone({
                    "points": coords,
                    "strokeWidth": 6,
                    "point_ids": [self.point2id[p1], self.point2id[p2]],
                    "link": link
                });

                link_group.add(hit_area, link)
                self.link_layer.add(link_group);
              }
           }
        }
      }
    }

    self.link_layer.draw();
    console.timeEnd("LINKS");
  }

  ChemSpace.prototype._get_links_feature = function(){
    var self = this;
    console.time("LINK COUNT [internal]");
    var p1, p2, links, link, link_list;
    var point2count = {};
    var done = {};
    self._links = false;

    for(var i = 0, len = self.points_len; i < len; i++){
      links = self.data.points[self.point_ids[i]].links;

      if(links !== undefined && links.length > 0){
        self._links = true;
        break;
      }
    }

    if(self._links){
        for(var i = 0, len = self.points_len; i < len; i++){
          p1 = self.point_ids[i];

          if(point2count[p1] === undefined){
            point2count[p1] = 0;
          }

          links = self.data.points[p1].links;

          if(links !== undefined && links.length > 0){

            for(var l = 0, len_links = links.length; l<len_links; l++){
              p2 = links[l][0];

              if(point2count[p2] === undefined){
                point2count[p2] = 0;
              }

              link = [p1,p2].sort().join("_");

              if(done[link] === undefined){
                done[link] = true;
                point2count[p1]++;
                point2count[p2]++;
              }
            }
          }
        }

        self.link_count = Object.keys(done).length;

        if(self.link_count > 0){
          for(var i = 0, len = self.points_len; i < len; i++){
            p1 = self.point_ids[i];
            self.data.points[p1][self.keys.features].push(point2count[p1]);
          }
          self.data.feature_names.push("Link count [internal]");
        }
    }
    else{
        self.settings.links.draw = false;
    }

    console.timeEnd("LINK COUNT [internal]");
  }

  ChemSpace.prototype._get_value_on_scale = function(value, prop){
    var self = this, result;
    var min_max_middle = self._prop2settings[prop].values;
    var scale = self.settings[prop].scale;

    if(value <= min_max_middle.min || value === null || value === undefined){
      return self.settings[prop].scale.min;
    }

    if(value >= min_max_middle.max){
      return self.settings[prop].scale.max;
    }
    if(value > min_max_middle.min && value <= min_max_middle.middle){
      radius = self._hack_round(scale.min + (value-min_max_middle.min)/(min_max_middle.middle-min_max_middle.min)*(scale.middle-scale.min));
    }
    else{
      radius = self._hack_round(scale.middle + (value-min_max_middle.middle)/(min_max_middle.max-min_max_middle.middle)*(scale.max-scale.middle));
    }
    return radius;

  }


  ChemSpace.prototype._draw_property_menu = function(property="point_size"){
      var self = this;
      if(!self.settings.navigation_toggle[property] || !self.settings.navigation_toggle.header|| (property == "link_width" && !self._links)){
        return;
      }

      if(self.target_element.find(".chemspacejs-navigation_div").length === 0){
        self.target_element.append($("<div class='chemspacejs-navigation_div'></div>")
            .css({
                "position": "absolute",
                "top": self.header_height - 70,
                "right": self.right_margin,
                "max-width": self.stage.getWidth() - self.left_margin - self.right_margin,
                "z-index": 10,
                "display": "grid",
                "grid-template-columns": "auto auto",
                "font-size": "0.8rem",
                "gap": "0.5rem"
            })
        );

      }
      var settings = self.settings[property];
      var font_size = 12;
      var margin = 5;
      var max_height = (self.settings.point_size.scale.max*2 > self.settings.link_width.scale.max)?self.settings.point_size.scale.max*2:self.settings.link_width.scale.max;
      if(max_height < 15){
        max_height = 15;
      }
      var order = ["min", "middle","max"];
      var space = 10;
      
      if((Number.isInteger(settings.index) && settings.index !== "category")||property == "link_width"){
        var menu_div = $("<div data-property='" + property + "'></div>")
            .css({"display": "flex", "flex-direction": "column", "align-items": "center", "padding": "5px 8px 0px 8px", "order": self._settings[property].order});
        var label = $("<div class='chemspacejs-label'>" + ((property === "link_width")?"Edge width":self.data.feature_names[self.settings[property].index]) +"</div>")
            .css({"font-weight": "bold"});
        menu_div.append(label);
        var values_div = $("<div></div>").css({"display": "flex", "justify-content": "space-between", "align-items": "center", "width": "100%", "font-size": "0.7rem"});
        
        if(property === "color"){
            var shape_div = $("<div class='chemspacejs-shape_div'></div>")
                .css({"display": "flex", "align-items": "center", "height": max_height, "margin-top": 3, "width": "100%"});
            var shape = self._settings[property].menu_shape.clone()
                .css({"background": "linear-gradient(to right, " + self._get_color_for_value(0, 0, 1, 0.5, self.settings.color.scale) + ", " + self._get_color_for_value(0.5, 0, 1, 0.5, self.settings.color.scale) + ", " + self._get_color_for_value(1, 0, 1, 0.5, self.settings.color.scale) + ")"});
            menu_div.append(shape_div.append(shape));

            for(var i=0, len=order.length; i<len; i++){
              var value_div = $("<div class='chemspacejs-value_div'></div>")
                .css({"display": "flex", "flex-direction": "column", "align-items": "center", "justify-content": "space-between", "padding": "3px"});

              var value = $("<div>" + ((settings.params[order[i]] === undefined||settings.value_type === "percentile")?self._prop2settings[property].values[order[i]]:settings.params[order[i]]) + "</div>")
               .css({"margin-top": 5});
              values_div.append(value_div.append(value));
            }
        }
        else{
            values_div.css({"margin-top": "3px"});

            for(var i=0, len=order.length; i<len; i++){
              var value_div = $("<div class='chemspacejs-value_div'></div>")
                .css({"display": "flex", "flex-direction": "column", "align-items": "center", "justify-content": "space-between", "padding": "3px"});
              var shape_div = $("<div class='chemspacejs-shape_div'></div>")
                .css({"display": "flex", "align-items": "center"});
              var shape = self._settings[property].menu_shape.clone();

              if(property == "point_size"){
                shape.css({"width": settings.scale[order[i]]*2, "height": settings.scale[order[i]]*2, "border-radius": settings.scale[order[i]]});
              }
              else{
                shape.css({"border-bottom-width": settings.scale[order[i]]});
              }
              value_div.append(shape_div.append(shape));
              var value = $("<div>" + ((settings.params[order[i]] === undefined||settings.value_type === "percentile")?self._prop2settings[property].values[order[i]]:settings.params[order[i]]) + "</div>")
               .css({"margin-top": 5});
              values_div.append(value_div.append(value));
            }
        }

        menu_div.append(values_div);
        self.target_element.find(".chemspacejs-navigation_div").append(menu_div);
        self.target_element.find(".chemspacejs-shape_div").css("height", max_height);

        menu_div.on("click", function(evt){self._property_menu_click(this, evt)});
        menu_div.on("mouseover", function(evt){self._icon_mouseover(evt, this)});
        menu_div.on("mouseout", function(evt){self._icon_mouseout(evt, this)});
        self._css_hover(menu_div);
      }
  }

  ChemSpace.prototype._property_menu_click = function(elm, evt){
    var self = this;
    var i, option, key, value;
    var property = $(elm).attr("data-property");
    var form_name = property+"_form";
    var order = ["min", "middle", "max"];
    var form = self.target_element.find("."+form_name);
    var overlay = self._draw_target_overlay();
    var settings = self.settings[property];

    if(form.length){
      form.fadeIn();
    }
    else{
      var form = self.html_ref.menu.clone().addClass(form_name)
        .css({"top": self.header_height, "right": self.right_margin, "max-width": 200});
      
      self.target_element.append(form);

      overlay.click(function(){
        form.fadeOut();
        overlay.fadeOut();
      });

      form.append(self.html_ref.menu_subtitle.clone().text(self._settings[property].property_label + " scale:"));

      var scale_div = self._min_max_middle_inputs(settings.scale, "scale", property).addClass("scale");
      form.append(scale_div);

      var value_type = $("<div></div>").css({"margin-bottom": 5});
      value_type.append(self.html_ref.menu_subtitle.clone().text(self._settings[property].property_label + " by:"));

      var value_type_select = $("<select name='value_type'>\
        <option value='percentile' " + ((settings.value_type == "percentile")?"selected":"") + ">Percentiles</option>\
        <option value='value' " + ((settings.value_type == "value")?"selected":"") + ">Values</option>\
        </select>").css({"padding": 3, "background-color": "white"});
      value_type.append(value_type_select);
      form.append(value_type);

      var values_div = self._min_max_middle_inputs((settings.value_type == "percentile")?settings.params:self._prop2settings[property].values, "params")
        .addClass("values");
      form.append(values_div);

      var redraw_button = self.html_ref.redraw_button.clone();
      form.append(redraw_button);
      form.find("input, select").css({"font-family": self.settings.font.fontFamily});
      self._css_hover(redraw_button);

      form.delegate(".redraw_button", "click", function(evt){
        var settings = {};
        var settings_fieldset = $(this).parents("."+form_name).find("input, select");

        settings_fieldset.each(function(){
            option = $(this);
            key = option.attr("name");
            value = option.val();
            datatype = option.attr("data-type");

            if(key == 'value_type'){
              self.settings[property][key] = value;
            }
            else if(value != ""){
              if(value[0] === "#"){
                var c = self._hex2rgb(value);
                var rgb = "rgb(" + c.r + "," + c.g + "," + c.b +")";
                self.settings[property][datatype][key] = rgb;
                self._color_dicts[rgb] = c;
              }
              else{
                self.settings[property][datatype][key] = parseFloat(value);
              }
            }
        });

        if(self.settings[property].value_type == "percentile"){
          self._percentile_alert(self.settings[property].params);
        }

        self.redraw_points();
        self._update_property_menu(property);
        self.navigation_layer.draw();
        overlay.trigger('click');
      });
    }
  }

  ChemSpace.prototype._update_property_menu = function(property){
    var self = this;
    var menu = self.target_element.find("[data-property='" + property + "']");

    // if(property === "color"){
    //  console.log(menu)
    // }

    if(menu.length !== 0){
      menu.remove();
    }
    self._draw_property_menu(property);
  }

  ChemSpace.prototype._draw_menu = function(){
    var self = this;
    self.navigation_layer = new Konva.Layer();

    if(self.settings.navigation_toggle.axis){
      self._draw_axis();
    }

    self._draw_axis_labels();
    self.stage.add(self.navigation_layer);
    var x = self.left_margin+40;
    var y = self.header_height-40;
    var scale = 1;
    if(!self.settings.navigation_toggle.header){
        x = self.left_margin + 5;
        y = self.header_height + 5;
        scale = 0.7;
    }
    
    if(self.settings.navigation_toggle.export_button){
      var export_icon = self._get_icon_group("export_icon", "Export\nin png format", {x: x, y: y});
      self.navigation_layer.add(export_icon);
      x = x+40;
    }

    self.refresh_icon = self._get_icon_group("refresh_icon", "Refresh", {x: x, y: y, scale: {x: scale, y: scale}});    
    self.navigation_layer.add(self.refresh_icon);
    self.refresh_icon.hide();

    if(!self.settings.navigation_toggle.header){
        return;
    }
    
    var navigation_selects = [
      {"name": "x", "label":"X axis"},
      {"name": "y", "label": "Y axis"},
      {"name": "color", "label": "Color"},
      {"name": "point_size", "label": "Point size"}
    ];

    self._draw_property_menu("color");
    self._draw_property_menu("point_size");
    self._draw_property_menu("link_width");

    var options_icon = self._get_icon_group("options_icon", "Options", {x: self.left_margin-3, y: y});
    self.navigation_layer.add(options_icon);


    if(self.target_element.find(".navigation").length == 0){
      var main_element = self.html_ref.menu.clone().addClass("navigation")
        .css({"top": 10, "left": 80, "display": "none"});

      var menu = $("<div></div>")
      main_element.append(menu);
      self.target_element.append(main_element);

      for(var i = 0, len=navigation_selects.length; i<len; i++){
        var current = navigation_selects[i];
        var label = self.html_ref.menu_subtitle.clone().text(current.label+":");
        var select = $("<select></select>").attr("data-name", current.name)
          .css({"padding": 3, "background-color": "white"});

        if(current.name === "point_size"){
          select.css("margin-bottom", 10);
        }

        for(var j = 0, len_2 = self.data.feature_names.length; j<len_2; j++){
          select.append($("<option value='" + j + "'>" + self.data.feature_names[j] + "</option>"));
        }

        if(current.name === "color" || current.name === "point_size"){
          select.append("<option value='category'>By category</option>");
          select.val((self.settings[current.name].index !== false)?self.settings[current.name].index:"category");
        }
        else{
          select.val(self.settings.coordinates[current.name]);
        }

        var div = $("<div></div>").css({"margin-bottom": 5});
        div.append(label, select);
        menu.append(div);
      }

      get_slider_val = function(value, slider){
        var target_id = $(slider).attr("id") + "_value";
        $("#" + target_id).text(value);
      }

      var slider_defs = [
        {"name": "resolution", "min": 0, "max": 50, "value": self.settings.resolution, "label": "Resolution", "type": "resolution", "step": 2},
      ];

      for(var i = 0, len=slider_defs.length; i<len; i++){
        var def = slider_defs[i];
        menu.append($("<div>\
          <label for=" + def.name + "_slider>" + def.label + ":</label>\
          <div id='" + self.settings.target + "_" + def.name + "_slider_value'>" + def.value + "</div>\
          </div>").css({"display": "flex", "justify-content": "space-between"}));

        menu.append($("<input oninput='get_slider_val(value, this)' type='range' step='5' class='slider' value="+ def.value +" min=" + def.min +" max=" + def.max + " data-type='" + def.type + "' id='" + self.settings.target + "_" + def.name + "_slider'/>")
          .css({"width": "100%", "padding": "5px 0px", "margin": 0}));
      }

      menu.append($("<div><input type='checkbox' name='links' checked/><div>Links</div></div>").css({"display": "flex", "align-items": "center", "margin-top": 5}));
      var name2layer = {
        "links": self.link_layer
      }

      menu.on("change", "input[type='checkbox']", function(){
        var state = $(this).prop("checked");
        var name = $(this).attr("name");
        var value = 1;
        self.settings[name].draw = state;

        if(!state){
          value = 0;
          $(this).attr("data-last_settings", [self.settings.coordinates.x, self.settings.coordinates.y].join("_"));
          $(name2layer[name].getCanvas()._canvas).hide();
        }
        else{
          $(name2layer[name].getCanvas()._canvas).show(); 
        }

        if(value == 1 && ([self.settings.coordinates.x, self.settings.coordinates.y].join("_") != $(this).attr("data-last_settings") || $(this).attr("data-last_x") !== self.selection_range.x.join("_") || $(this).attr("data-last_y") !== self.selection_range.y.join("_"))){
          $(this).attr("data-last_x", self.selection_range.x.join("_"));
          $(this).attr("data-last_y", self.selection_range.y.join("_"));
          self._destroy_children([self.link_layer]);
          self._draw_links();
        }

        
        self._sort_layers();
      });

      var redraw_button = self.html_ref.redraw_button.clone();
      self._css_hover(redraw_button);
      menu.append(redraw_button);
      menu.find("input, select").css({"font-family": self.settings.font.fontFamily});

      main_element.find("label").css({"font-weight": "bold"});
      main_element.find("input[type='checkbox']").css({"margin-right": 5});

      $(".slider").on("change", function(){
        var slider_type = $(this).attr("data-type");
        var value = parseInt($(this).val());
        self.settings.resolution = value;
        self._prop2settings.coordinates.x = -1;
      });

      menu.on("change", "select", function(){
        var name = $(this).attr("data-name");
        var value = $(this).val();

        if(name === 'color' || name === 'point_size'){
          self.settings[name].index = (isNaN(parseInt(value)))?value:parseInt(value);
          console.log(self.settings[name].index)
        }
        else{
          self.settings.coordinates[name] = (isNaN(parseInt(value)))?value:parseInt(value);
        }

        if(name === "color"){
            if(self.settings.color.value_type == "value"){
              self.settings.color.params = {};
              var color_settings = self.target_element.find(".color_settings");

              color_settings.find(".max").val(self._prop2settings.color.values.abs_max);
              color_settings.find(".min").val(self._prop2settings.color.values.abs_min);
              color_settings.find(".middle").val(self._prop2settings.color.values.middle);
            }
        }
      });

      redraw_button.on("click", function(){
        self.redraw_points();
        self._draw_legend();
        self._draw_axis_labels();
        self._update_property_menu("color");
        self._update_property_menu("point_size");
        self._update_property_menu("link_width");
        self.navigation_layer.draw();
        self._draw_compounds();
        self.target_element.find(".target_overlay").trigger('click');
      });

    }


  }

  ChemSpace.prototype._draw_axis = function(){
    var self = this;
    var x = self.selection_overlay_rect.attrs.x;
    var y = self.selection_overlay_rect.attrs.y;
    var top = self.objects_ref.border_line.clone({points:[x, y, x + self.selection_overlay_rect.attrs.width, y]});
    var right = self.objects_ref.border_line.clone({points:[x + self.selection_overlay_rect.attrs.width, y, x + self.selection_overlay_rect.attrs.width, y + self.selection_overlay_rect.attrs.height]});
    var bottom = self.objects_ref.border_line.clone({points:[x, y + self.selection_overlay_rect.attrs.height, x + self.selection_overlay_rect.attrs.width, y + self.selection_overlay_rect.attrs.height]});
    var left = self.objects_ref.border_line.clone({points:[x, y + self.selection_overlay_rect.attrs.height, x, y]});
    self.navigation_layer.add(top, right, bottom, left);
  }

  ChemSpace.prototype._draw_axis_labels = function(){
    var self = this;

    if(self.settings.navigation_toggle.axis_labels == false){
      return;
    }

    if(self.x_axis_label_group !== undefined){
      self.x_axis_label_group.destroy();
      self.y_axis_label_group.destroy();
    }

    self.x_axis_label_group = new Konva.Group({listening: false});
    self.y_axis_label_group = self.x_axis_label_group.clone();

    var y = self.settings.height - self.footer_height + 5;
    var x_text, y_text, x_min, x_max, y_min, y_max;

    x_min = self.objects_ref.navigation_text.clone({
      x: self.left_margin + 5,
      y: y,
      text: self._hack_round_float(self.coordinate_ranges["x"][0], 2).toString(),
    });

    x_label = self.objects_ref.navigation_text.clone({
      y: y,
      text: self.data.feature_names[self.settings.coordinates.x].toString(),
      fontStyle: "italic"
    });

    x_label.setAttrs({
      x: (self.settings.width - self.left_margin - self.right_margin)/2 + self.left_margin - x_label.width()/2,
    });

    x_max = self.objects_ref.navigation_text.clone({
      y: y,
      text: self._hack_round_float(self.coordinate_ranges["x"][1], 2).toString(),
    });

    x_max.setAttrs({
      x: self.settings.width - self.right_margin - x_max.width() - 5,
    });
    self.x_axis_label_group.add(x_label, x_min, x_max);

    var x2 = self.left_margin - self.settings.font.fontSize - 5;
    y_min = self.objects_ref.navigation_text.clone({
      x: x2,
      y: self.settings.height - self.footer_height - 10,
      text: self._hack_round_float(self.coordinate_ranges["y"][0], 2).toString(),
      rotation: -90,
    });

    y_label = self.objects_ref.navigation_text.clone({
      x: x2,
      text: self.data.feature_names[self.settings.coordinates.y].toString(),
      rotation: -90,
      fontStyle: "italic"
    });

    y_label.setAttrs({
      y: (self.settings.height - self.footer_height - self.header_height)/2 + self.header_height + y_label.width()/2,
    });

    y_max = self.objects_ref.navigation_text.clone({
      x: x2,
      text: self._hack_round_float(self.coordinate_ranges["y"][1], 2).toString(),
      rotation: -90,
    });

    y_max.setAttrs({
      y: self.header_height + y_max.width() + 5,
    });
    self.y_axis_label_group.add(y_label, y_min, y_max);
    self.navigation_layer.add(self.x_axis_label_group, self.y_axis_label_group);
  };

  ChemSpace.prototype._feature_switch = function(){
    var self = this;
    var feature_element = self.target_element.find(".navigation");
    if(feature_element.is(":visible")){
      feature_element.fadeOut();
    }
    else{
      var overlay = self._draw_target_overlay();
      feature_element.fadeIn();

      overlay.click(function(){
        feature_element.fadeOut();
        overlay.fadeOut();
      });
    }
  }

  ChemSpace.prototype._draw_legend = function(){
    var self = this;
    var shapes_css = {
        "circle": {
            "min-width": 16, "min-height": 16, "max-width": 16, "max-height": 16, "border-radius": 8
        },
        "square": {
            "min-width": 14, "min-height": 14, "max-width": 14, "max-height": 14
        },
        "triangle": {
            "max-width": 0, "max-height": 0, "border-left": "8px solid transparent", "border-right": "8px solid transparent", "border-bottom": "16px solid red"
        },
        "rhombus": {
            "min-width": 12, "min-height": 12, "max-width": 12, "max-height": 12, "transform": "rotate(45deg)"
        },
        "triangle-down": {
            "max-width": 0, "max-height": 0, "border-left": "8px solid transparent", "border-right": "8px solid transparent", "border-bottom": "16px solid red", "transform": "rotate(180deg)"
        },
    };

    if(self.settings.navigation_toggle.header && self.settings.navigation_toggle.point_count){
        if(self.point_count !== undefined){
          self.point_count.destroy();
        }

        self.point_count = self.objects_ref.legend_text.clone({
          "text": Object.keys(self.point2coord).length + " / " + self.points_len,
          "fontSize": self.settings.font.fontSize,
          "y": self.header_height + 7,
        });

        var x = self.settings.width - self.right_margin - self.point_count.getWidth() - 7;
        self.point_count.setAttr("x", x);
        self.navigation_layer.add(self.point_count);
    }

    if(self.categories.length < 2){
      return;
    }

    var legend_div = self.target_element.find(".chemspacejs-legend_div");
    var first = false;
    var ls = self.settings.categories.legend;
    if(legend_div.length){
        var legend_elm = legend_div.find(".chemspacejs-legend");
        legend_elm.html("");
    }
    else{
        legend_div = $("<div class='chemspacejs-legend_div'></div>").css({"position": "absolute", "top": 20, "overflow": "auto"});
        var legend_elm = $("<div class='chemspacejs-legend'></div>").css({"font-weight": ls["fontWeight"], "font-size": ls["fontSize"], "font-family": ls["fontFamily"], "color": ls["fill"]});
        legend_div.css({"right": 0, "max-width": ls["width"], "width": ls["width"], "min-width": ls["width"], "max-height": self.settings.height});
    }

    for(var i = 0, len=self.categories.length; i < len; i++){
      category = self.categories[i];
      var l = $("<div class='chemspacejs-category' data-ci='" + i + "'></div>").css({"display": "flex", "align-items": "center"});
      if(!category.status){
        l.css("opacity", 0.3);
      }
      var color = (self.settings.color.index == "category")?category.color:self.settings.shapes.circle.fill;
      var fill_color = (category.status)?color:self.settings.shapes.circle.fill;
      var shape = (category.shape !== undefined)?category.shape:"circle";
      var shape_div = $("<div class='chemspacejs-legend_shape_div'></div>").css({"min-width": 20, "display": "flex", "justify-content": "center"});
      var shape_elm = $("<div class='chemspacejs-legend_shape'></div>").css((shapes_css[shape] !== undefined)?shapes_css[shape]:shapes_css.circle);
      shape_div.append(shape_elm);

      if(shape === "triangle" || shape === "triangle-down"){
        shape_elm.css({"border-bottom-color": color});
      }
      else{
        shape_elm.css({"background-color": fill_color, "border-color": color});
      }

      l.append(shape_div);
      var label = $("<div>" + category.label + " (" + self.category2count[i] + "/" + category.points.length + ")" + "</div>").css({"padding": "3px 8px", "font-family": ls["fontFamily"]});
      l.append(label);
      legend_elm.append(l);
    }
    legend_div.append(legend_elm);
    self.target_element.append(legend_div);

    self._css_hover(legend_elm.find(".chemspacejs-category"));

    legend_elm.find(".chemspacejs-category").on("click", function(evt){
        self._category_legend_click($(this), evt);
    });

    legend_elm.find(".chemspacejs-category").on("mouseover", function(evt){
        self._category_legend_mouseover($(this), evt);
    });

    legend_elm.find(".chemspacejs-category").on("mouseout", function(evt){
        self._category_legend_mouseout($(this), evt);
    });
  }

  ChemSpace.prototype._process_categories = function(categories){
    var self = this;
    self.category2layer = {};
    self.categories = [];
    var shape_index = 0;

    if(categories !== undefined && self.settings.categories.draw){
        if(self.settings.categories.order !== false){
          for(var i = 0, len = categories.length; i<len; i++){
            var c = categories[i];
            var ci = self.settings.categories.order.indexOf(c.label);
            if(ci !== -1){
              categories[i].order = ci;
            }
          }
        }
      
      categories = categories.sort(function(a, b){return a.order > b.order});
      for(var i = 0, len = categories.length; i < len; i++){
        var c = categories[i];
        self.category2layer[i] = new Konva.Layer({"category_index": i, "id": self.settings.target + "_category_" + i, "label": c.label});
        if(c.color === undefined){
          c.color = self.settings.colors[i];
        }

        if(c["shape"] === undefined){
          if(shape_index == self.settings.shapes.order.length){
            shape_index = 0;
          }
          c.shape = self.settings.shapes.order[shape_index];
          shape_index++;
        }

        if(c["radius"] === undefined){
          c.radius = self.settings.shapes.circle.radius;
        }

        if(c.points === undefined && c.objects !== undefined){
          c.points = self._translate_objects_to_points(c.objects);
        }

        c.status = true;
        self.categories.push(c);
      }
    }
  }

  ChemSpace.prototype.get_points_for_object = function(object_id){
    var self = this;
    return self._translate_objects_to_points([object_id]);
  }

  ChemSpace.prototype.get_data_for_point = function(point_id){
    var self = this;
    var data = self.data.points[point_id];
    data[self.keys.smiles] = (self.data.compounds[point_id] !== undefined)?self.data.compounds[point_id][self.keys.smiles]:null;
    data.color = self._id2object[self.point2id[point_id]].attrs.fill;
    return data;
  }

  ChemSpace.prototype._calculate_coordinates = function(point_ids){
    var self = this;
    console.log(self.data.feature_names)

    if(self.use_cache && point_ids == undefined && [self._prop2settings.coordinates.x, self._prop2settings.coordinates.y].join("_") == [self.settings.coordinates.x, self.settings.coordinates.y].join("_")){
      console.log("COORDINATES CACHE");
      return true;
    }
    console.time("CALCULATING COORDINATES");
    self._prop2settings.coordinates = $.extend(true, {}, self.settings.coordinates);
    self.point2coord = {};
    self.id2data = {};
    var coord, key, data;

    if(point_ids === undefined){
      var point_ids = self.point_ids;
      data = self.data.points;
    }
    else{
      data = {};
      for(var i = 0, len = point_ids.length; i<len; i++){
        key = point_ids[i];
        data[key] = self.data.points[key];
      }
    }

    var coord, key, item, xs=[], ys=[], x_point_index, y_point_index;
    var x_index = self.settings.coordinates.x;
    var y_index = self.settings.coordinates.y;
    
    for(var i = 0, len=point_ids.length; i<len; i++){
        item_data = data[point_ids[i]][self.keys.features];
        xs.push(item_data[x_index]);
        ys.push(item_data[y_index]);
    };

    self.coordinate_ranges = {"x": self._get_min_max(xs), "y": self._get_min_max(ys)};

    var x_border = (self.coordinate_ranges["x"][1] - self.coordinate_ranges["x"][0])*0.01;
    var max_x = self.coordinate_ranges["x"][1]+x_border;
    var min_x = self.coordinate_ranges["x"][0]-x_border;

    var y_border = (self.coordinate_ranges["y"][1] - self.coordinate_ranges["y"][0])*0.01;
    var max_y = self.coordinate_ranges["y"][1]+y_border;
    var min_y = self.coordinate_ranges["y"][0]-y_border;

    if(self.settings.compounds.draw){
      self.compound_size = self.settings.compounds.size;
      self.atom_size = self.settings.compounds.atom_size + (1 - point_ids.length/self.points_len)*(self.settings.compounds.atom_size - self.settings.compounds.atom_size);
    }

    var margin = ((self.settings.compounds.draw && point_ids.length <= self.settings.compounds.limit)?self.compound_size/2:(self.settings.point_size.scale.max/2)) + 15;
    var origin_x = self._hack_round(self.left_margin + margin);
    var origin_y = self._hack_round(self.header_height + margin);


    var width = self.settings.width - self.left_margin - self.right_margin - 2*margin;
    var height = self.settings.height - self.header_height - self.footer_height - 2*margin;
    var x_coord, y_coord;
    var res = (self.settings.resolution == 0)?1:self.settings.resolution;
    var x_part = self._hack_round(width/parseInt(width/res));
    var y_part = self._hack_round(height/parseInt(height/res));
    var x_shift = parseInt(x_part/2);
    var y_shift = parseInt(y_part/2);

    self.point_index = {};
    
      for(var i = 0, len = point_ids.length; i < len; i++){
        key = point_ids[i];
        coord = data[key][self.keys.features];

        x_coord = width*(coord[x_index]-min_x)/(max_x-min_x);
        x_point_index = origin_x + (x_coord/x_part << 0)*x_part+x_shift;

        y_coord = height*(coord[y_index]-min_y)/(max_y-min_y);
        y_point_index = origin_y + height - (y_coord/y_part << 0)*y_part - y_shift;

        self.point2coord[key] = [x_point_index, y_point_index];

        if(self.point_index[x_point_index] === undefined){
          self.point_index[x_point_index] = {};
        }

        if(self.point_index[x_point_index][y_point_index] === undefined){
          self.point_index[x_point_index][y_point_index] = [key];
        }
        else{
          self.point_index[x_point_index][y_point_index].push(key);
        }
      }
    
    if(self.data.feature_names.indexOf("Point count [internal]") === -1){
      self.data.feature_names.push("Point count [internal]");
    }
    
    let point_count_index = self.data.feature_names.indexOf("Point count [internal]");

    let x_keys = Object.keys(self.point_index);
    for(let i = 0, x_len=x_keys.length; i < x_len; i++){
        let x_key = x_keys[i];
        let y_keys = Object.keys(self.point_index[x_key]);

        for(let j = 0, y_len=y_keys.length; j < y_len; j++){
          let y_key = y_keys[j];
          let point_ids = self.point_index[x_key][y_key];
          let point_count = point_ids.length;

          for(let p = 0, p_len=point_count; p<p_len; p++){
            let point_id = point_ids[p];
            self.data.points[point_id][self.keys.features][point_count_index] = point_count;
          }
        }
    }
    console.timeEnd("CALCULATING COORDINATES");;
    return false;
  }

  /*ChemSpace.prototype._draw_paths = function(){
    var self = this;
    var path;
    console.time("PATHS");

    for(var i = 0, len = self.data.paths.length; i < len; i++){
      var points = [];
      for(var p = 0, len_points = self.data.paths[i]["points"].length; p<len_points; p++){
        var point = self.point2coord[self.data.paths[i]["points"][p]];
        if(point !== undefined){
          points = points.concat(point);
        }
      }
      var color = self.data.paths[i].color;
      path = self.objects_ref.path.clone({"points": points, "id": self.data.paths[i]["label"], "path_index": i, "stroke": (color !== undefined)?color:self.settings.path.stroke});
      self.path_layer.add(path);
      // path.moveToBottom();
    }
    self.path_layer.draw();
    console.timeEnd("PATHS");
  }*/

  ChemSpace.prototype._get_id2categories = function(){
    var self = this;
    
    var id2categories = {}, key;
    var create_other_category = false;

    if(self.settings.categories.draw){
        for(var i = 0, len=self.points_len; i<len; i++){
          id2categories[self.point_ids[i]] = [];
        }

        for(var i = 0, len=self.categories.length; i < len; i++){
          var category = self.categories[i];

          for(var j = 0, len_2=category.points.length; j < len_2; j++){
            if(id2categories[category.points[j]] !== undefined){
              id2categories[category.points[j]].push(i);
            }
          }
        }

        for(var i = 0, len=self.points_len; i<len; i++){
          if(id2categories[self.point_ids[i]].length == 0){
            create_other_category = true;
            break;
          }
        }

        if(create_other_category){
          self.category2layer[self.categories.length] = new Konva.Layer();
          self.categories.push({"status": true, "label": "other", "color": self.settings.shapes.circle.fill, "radius": self.settings.shapes.circle.radius, "points": []});
          var other_index = self.categories.length - 1;

          for(var i = 0, len=self.points_len; i<len; i++){
            key = self.point_ids[i];

            if(id2categories[key].length == 0){
              self.categories[other_index].points.push(key);
              id2categories[key].push(other_index);
            }
          }
        }
    }
    else{
        self.category2layer[0] = new Konva.Layer();
        self.categories.push({"status": true, "label": "other", "color": self.settings.shapes.circle.fill, "radius": self.settings.shapes.circle.radius, "points": self.point_ids});
        for(var i = 0, len=self.points_len; i<len; i++){
          id2categories[self.point_ids[i]] = [0]
        }
    }

    return id2categories;
  }

  /************************************
    Event functions
  ************************************/

  ChemSpace.prototype._icon_click = function(evt){
    var self = this;
    var icon_name = evt.target.attrs.icon_name;

    if(icon_name === "zoom_icon"){
      self._zoom_icon_click();
    }
    else if(icon_name === "options_icon"){
      self._feature_switch();
    }
    else if(icon_name === "export_icon"){
      self.export_image();
    }
    else if(icon_name === "refresh_icon"){
      self.refresh_icon.hide();
      self.redraw_points(self.point_ids);
      self._draw_legend();
      self._draw_axis_labels();
      self.navigation_layer.draw();
      self._draw_compounds();
    }
    else if(icon_name === "color_scale_icon"){
      self._color_scale_click();
    }
    /*else if(icon_name === "point_size"||icon_name === "link_width"){
      self._property_menu_click(evt);
    }*/
  }

  ChemSpace.prototype.export_image = function(action, filename){
    var self = this;
    self.hover_layer.destroyChildren();
    self.hover_layer.draw();

    html2canvas(self.target_element[0], {backgroundColor: null}).then(canvas => {
        // document.body.appendChild(canvas);
        var a = document.createElement('a');
        a.href = canvas.toDataURL("image/png").replace("image/png", "image/octet-stream");
        a.download = 'chemspacejs.png';
        a.click();
    });
  }

  ChemSpace.prototype._zoom_icon_click = function(){
    var self = this;
    self.hover_layer.destroyChildren();
    self.hover_layer.draw();
    self.selection_rect.destroy();
    self.selection_rect = false;
    self.zoom_icon.destroy();
    self.selection_layer.draw();
    if(self.refresh_icon !== undefined){
        self.refresh_icon.show();
    }

    self.redraw_points(self.current_selection);
    self._draw_legend();
    self.navigation_layer.draw();
    self._draw_axis_labels();
    self.navigation_layer.draw();
    self.navigation_layer.moveToTop();
    self._draw_compounds();
  }

  ChemSpace.prototype._point_click = function(evt) {
    var self = this;
    var point = evt.target;
    var point_ids = point.attrs.point_ids;
    var object_ids = [];

    for(var i = 0, len = point_ids.length; i<len; i++){
      object_ids = object_ids.concat(self.data.points[point_ids[i]][self.keys.object_ids]);
    }
    
    self.events.point_click(object_ids, evt);
  }

  ChemSpace.prototype._get_point_tooltip = function(evt){
    var self = this;
    var select, row;
    var obj_class = evt.target.className;
    
    var tooltip_wrapper = $("<div class='point_tooltip_wrapper'></div>")
      .css({
        "position": "fixed",
        "display": "none",
        "z-index": 100,
        "align-items": "center"
      });

    if(obj_class === "Line"){
        var point_ids = evt.target.attrs.point_ids;
        var tooltip = self.tooltip_templates.default(point_ids, tooltip_wrapper);
    }
    else{
        var point_ids = [evt.target.attrs.id];
        var tooltip = self.tooltip_templates.default(point_ids, tooltip_wrapper);
    }    

    return tooltip;
  }

  ChemSpace.prototype._get_points_data = function(point_ids){
    var self = this;

    if(point_ids.length == 1){
      var point_data = self.data.points[point_ids[0]][self.keys.features];
    }
    else{
      var point_data = [];
      var pids_len = point_ids.length;

      for(var i = 0; i<self.data.feature_names.length; i++){
        var values = [];

        for(var j = 0; j<pids_len; j++){
          values.push(self.data.points[point_ids[j]][self.keys.features][i]);
        }
        values.sort(self._sort_number_ascending);
        point_data.push((values[0] !== null)?values[0] + " - " + values[pids_len-1]: "-");
      }
    }
    return point_data;
  }

  ChemSpace.prototype._point_mouseover = function(evt) {
    var self = this;

    var tooltip_check = self.target_element.find(".point_tooltip");

    if(tooltip_check.length){
      tooltip_check.remove();
    }

    var point_ids = self.id2points[evt.target.attrs.id];
    var coords = self.point2coord[point_ids[0]];
    var color = evt.target.fill();
    var tooltip = self.events.point_tooltip(self.point_index[coords[0]][coords[1]], color, evt);
    
    self.target_element.append(tooltip);
    self._place_tooltip(tooltip, coords, evt);

    if(self.settings.compounds.draw){
      self._draw_compound_in_tooltip(point_ids);
    }
    
    if(evt.target.className === 'Circle'||evt.target.className === 'RegularPolygon'){
      var hover_obj = evt.target.clone({"listening": false, "radius": evt.target.attrs.radius + 3});
    }
    else if(evt.target.className === 'Rect'){
     var hover_obj = evt.target.clone({
        "listening": false,
        "width": evt.target.attrs.width + 6,
        "height": evt.target.attrs.width + 6,
        "x": evt.target.attrs.x - 3,
        "y": evt.target.attrs.y - 3,
      });
    }

    self.hover_layer.add(hover_obj);
    self.hover_layer.draw()
  }

  ChemSpace.prototype._place_tooltip = function(tooltip, coords, evt){
    var self = this;
    tooltip.css({
      "display": "flex"
    });
    var current_x = evt.evt.clientX;
    var current_y = evt.evt.clientY;
    var tooltip_height = tooltip.height();
    var x = (coords[0] > self.settings.width/2)?(current_x - tooltip.width() - 20):(current_x + 20);
    var y = current_y - tooltip_height/2;
    var target_top = self.target_element[0].getBoundingClientRect().top;

    if(y < target_top){
      y = target_top;
    }

    if((y + tooltip_height) > (target_top + self.settings.height)){
      y = target_top + self.settings.height - tooltip_height;
    }

    tooltip.css({
      "top": y,
      "left": x,
    });
  }

  ChemSpace.prototype._point_mouseout = function(evt) {
    var self = this;
    self.target_element.find(".point_tooltip_wrapper").remove();
    self.hover_layer.destroyChildren();
    self.hover_layer.draw();
  }

  ChemSpace.prototype._path_mouseout = function(evt) {
    var self = this;
    self.hover_layer.destroyChildren();
    self.hover_layer.draw();
  }

  ChemSpace.prototype._path_mouseover = function(evt) {
    var self = this;
    var path = self.stage.find("#" + evt.target.attrs.id)[0];
    self.hover_layer.add(path.clone({"opacity": 1}))
    self.hover_layer.draw()
  }

  ChemSpace.prototype._path_click = function(evt) {
    var self = this;
    var path = self.data.paths[evt.target.attrs.path_index];
    self.events.path_click(path, evt);
  }

  ChemSpace.prototype._link_mouseout = function(evt) {
    var self = this;
    self.target_element.find(".point_tooltip_wrapper").remove();
    self.hover_layer.destroyChildren();
    self.hover_layer.draw();
  }

  ChemSpace.prototype._link_mouseover = function(evt) {
    var self = this;
    var point_ids = evt.target.attrs.point_ids;
    var link = evt.target.attrs.link.clone({"listening": false, "strokeWidth": evt.target.attrs.link.strokeWidth() + 2});
    var weight = evt.target.attrs.weight;
    self.hover_layer.add(link);
    
    var p1 = self.stage.find("#"+point_ids[0])[0];
    if(p1 !== undefined){
        var point_1 = p1.clone({"listening": false, "radius": p1.attrs.radius + 2});
        self.hover_layer.add(point_1);
    }

    var p2 = self.stage.find("#"+point_ids[1])[0];
    if(p2 !== undefined){
        var point_2 = p2.clone({"listening": false, "radius": p2.attrs.radius + 2});
        self.hover_layer.add(point_2);
    }

    self.hover_layer.draw()

    var coords_1 = self.point2coord[self.id2points[point_ids[0]][0]];
    var coords_2 = self.point2coord[self.id2points[point_ids[1]][0]];
    var coords = [(coords_1[0]+coords_2[0])/2, (coords_1[1]+coords_2[1])/2];

    var color = evt.target.fill();
    var tooltip = self.events.point_tooltip(point_ids, color, evt);

    self.target_element.append(tooltip);
    self._place_tooltip(tooltip, coords, evt);

    if(weight != undefined){        
        var weight_div = $("<div><div>Weight</div></div>")
            .css({
                "position": "absolute",
                "display": "flex",
                "flex-direction": "column",
                "color": "white",
                "align-items": "center",
                "background-color": "gray",
                "padding": "5px 10px",
                "font-size": "90%",
            });
        weight_div.append($("<div>" + weight +"</div>").css({"font-weight": "bold"}));      
        tooltip.append(weight_div);
        weight_div.css("right", tooltip.width()/2 - weight_div.width() - 10);
        weight_div.css("top", -weight_div.height()/2 - 5);
    }

    if(self.settings.compounds.draw){
      self._draw_compound_in_tooltip([self.id2points[point_ids[0]][0], self.id2points[point_ids[1]][0]]);
    }
  }

  ChemSpace.prototype._link_click = function(evt) {
    var self = this;
    var point_ids = evt.target.attrs.point_ids;
    self.events.link_click([self.id2points[point_ids[0]], self.id2points[point_ids[1]]], evt);
  }

  ChemSpace.prototype._sort_layers = function(){
    var self = this;
    console.time("SORT LAYERS");
    self.main_layer.moveToTop();
    if(self.path_layer !== undefined){
      self.path_layer.moveToTop();
    }

    if(self.link_layer !== undefined && self.settings.links.draw){
      self.link_layer.moveToTop();
    }

    if(self.compounds_layer !== undefined && self.settings.compounds.draw){
      self.compounds_layer.moveToTop();
    }

    if(self.categories[self.categories.length-1].label == "other" && self.categories[self.categories.length-1].status){
      self.category2layer[self.categories.length-1].moveToTop();
    }
    for(var i = self.categories.length-1; i>=0; i--){
      if(self.categories[i].status){
        if(self.categories[i].label != "other"){
          self.category2layer[i].moveToTop();
        }
      }
    }

    if(self.highlight_layer !== undefined){
      self.highlight_layer.moveToTop();
    }

    self.hover_layer.moveToTop();
    self.selection_layer.moveToTop();
    console.timeEnd("SORT LAYERS")
  }

  ChemSpace.prototype._category_legend_click = function(category_elm, evt) {
    var self = this;
    var ci = category_elm.attr("data-ci");
    var category = self.categories[ci];
    var shape_elm = category_elm.find(".chemspacejs-legend_shape");
    var on_categories = [];
    if(category.status){
      category_elm.css("opacity", 0.3);
      category.status = false;
    }
    else{
      on_categories.push(category);
      category_elm.css("opacity", 1);
      category.status = true;
    }
    
    if(self.compounds_layer.getChildren().length > 0){
        self._update_visible_compounds();
    }

    self.events.category_legend_click(category, evt);
    self.navigation_layer.draw();
    self._sort_layers();
  }

  ChemSpace.prototype._update_visible_compounds = function(on_categories){
    var self = this;

    if(self.compounds_layer.getChildren().length == 0){
        return;
    }

    if(on_categories === undefined){
        var on_categories = self.categories.map(function(c, i){if(c.status){return i}});
    }
    var current_point_ids = Object.keys(self.point2coord);

    for(var i = 0, len=current_point_ids.length; i<len; i++){
        var opacity = self.settings.categories.off_opacity;
        var point_id = current_point_ids[i];
        for(var i1 = 0, len1=self.point_id2categories[point_id].length; i1<len1; i1++){
            if(on_categories.indexOf(self.point_id2categories[point_id][i1]) > -1){
                opacity = 1;
                break;
            }
        }
        self.point2img[point_id].setAttr("opacity", opacity);
    }
    self.compounds_layer.draw();
  }

  ChemSpace.prototype._category_legend_mouseover = function(category_elm, evt) {
    var self = this;
    var ci = parseInt(category_elm.attr("data-ci"));
    
    for(var i = 0, len=self.categories.length; i<len; i++){
        if(i === ci){
            $(self.category2layer[i].canvas.context.canvas._canvas).css("opacity", 1);
            self.category2layer[i].setAttr("listening", true).moveToTop();

        }
        else{
            $(self.category2layer[i].canvas.context.canvas._canvas).css("opacity", self.settings.categories.off_opacity);
            self.category2layer[i].setAttr("listening", false);
        }
    }
    self.hover_layer.moveToTop();
    self._update_visible_compounds([ci]);    
  }

  ChemSpace.prototype._category_legend_mouseout = function(category_elm, evt) {
    var self = this;
    self.navigation_layer.draw();

    for(var i = 0, len=self.categories.length; i<len; i++){
        if(self.categories[i].status){
            $(self.category2layer[i].canvas.context.canvas._canvas).css("opacity", 1);
            self.category2layer[i].setAttr("listening", true);
        }
        else{
            $(self.category2layer[i].canvas.context.canvas._canvas).css("opacity", self.settings.categories.off_opacity);
            self.category2layer[i].setAttr("listening", false);
        }
    }
    self._update_visible_compounds();
  }

  ChemSpace.prototype._icon_mouseover = function(evt, elm){
    var self = this;
    
    if(elm !== undefined){
        var target = $(elm);
        var label = self._settings[target.attr("data-property")].label;
        var y = self.header_height;
        var t_bcr = self.target_element[0].getBoundingClientRect();
        var e_bcr = elm.getBoundingClientRect();        
        var x = e_bcr.x - t_bcr.x;      
    }
    else{
        var target = evt.target;
        var label = target.getAttr("label");
        var x = target.getAttr("x");
        var y = target.getAttr("y");
        var height = target.getHeight();
        y = y+height;
    }

    var icon_tooltip = self.objects_ref.tooltip_label.clone({x: x,
        y: y
    });

    icon_tooltip.add(self.objects_ref.tooltip_tag.clone());
    icon_tooltip.add(self.objects_ref.tooltip_text.clone({text: label}));
    self.hover_layer.moveToTop();
    self.hover_layer.add(icon_tooltip);
    self.hover_layer.draw();
  }

  ChemSpace.prototype._icon_mouseout = function(evt){
      var self = this;
      self.hover_layer.destroyChildren();
      self.hover_layer.draw();
  }

  ChemSpace.prototype._get_min_max_middle = function(data, percentiles){
    var self = this;
    var i, len;
    var min_max_middle = {};

    if(percentiles == undefined){
      var percentiles = {"min": 0, "max": 100, "middle": 50};
    }

    var len = data.length - 1;
    data.sort(self._sort_number_ascending);

    min_max_middle["min"] = (percentiles.min > 0)?data[self._hack_round((len-1)*percentiles.min/100)]:data[0];
    min_max_middle["max"] = (percentiles.max < 100)?data[self._hack_round((len-1)*percentiles.max/100)]:data[len];;
    min_max_middle["middle"] = (percentiles.middle != 50)?data[self._hack_round((len-1)*percentiles.middle/100)]:data[self._hack_round((len-1)/2)];

    min_max_middle["median"] = data[self._hack_round((len-1)/2)];
    min_max_middle["abs_min"] = data[0];
    min_max_middle["abs_max"] = data[len];

    return min_max_middle;
  }

  ChemSpace.prototype._get_prop_settings = function(prop){
    var self = this;
    var i, len;
    var settings = self.settings[prop];
    var cs = self._prop2settings[prop];
    var current_scale = (typeof(cs.scale) != "object")?cs.scale:[cs.scale.min, cs.scale.middle, cs.scale.max, cs.resolution].join("_");
    var new_scale = (typeof(settings.scale) != "object")?settings.scale:[settings.scale.min, settings.scale.middle, settings.scale.max, settings.resolution].join("_");

    if(!self.use_cache || [current_scale, cs.value_type, cs.index, cs.params.min, cs.params.middle, cs.params.max].join("_") != [new_scale, settings.value_type, settings.index, settings.params.min, settings.params.middle, settings.params.max].join("_")){
      var prop2null = {
        "color": "#D2D2D2",
        "point_size": self.settings.point_size.scale.min,
        "link_width": 0.5
      };

      console.log(prop, "CALCULATION");
      self._prop2settings[prop] = $.extend(true, {}, settings);
      self._prop2settings[prop].resolution = self.settings.resolution;
      var v2v = {};

      if(self._prop2settings[prop].index == 'category'){
        for(var i = 0, len=self.categories.length; i<len; i++){
          v2v[i] = self._prop2fnc[prop].category(self.categories[i]);
        }
      }
      else{
        var keys = self.point_ids;
        var data = [];

        if(prop === "link_width"){
            for(var i = 0, len=keys.length; i<len; i++){
              var point = self.data.points[keys[i]];

              if(point.links !== undefined){
                for(var v = 0, len_v=point.links.length; v < len_v; v++){
                    var value = point.links[v][1]
                    data.push(value);
                    v2v[value] = null;
                }   
              }
            }
        }
        else{
            var prop_index = self.settings[prop].index;
            for(var i = 0, len=keys.length; i<len; i++){
              var value = self.data.points[keys[i]][self.keys.features][prop_index];
              // console.log(self.data.points[keys[i]])
              if(value !== null && value !== undefined){
                data.push(value);
                v2v[value] = null;
              }
            }
        }


        var len = data.length;
        data.sort(self._sort_number_ascending);
        var min_max_middle = self._get_min_max_middle(data, self.settings[prop].params);
        // console.log(min_max_middle, data)
        if(self.settings[prop].value_type === "value"){
          min_max_middle["min"] = (settings.params.min !== undefined)?settings.params.min:min_max_middle["abs_min"];
          min_max_middle["max"] = (settings.params.max !== undefined)?settings.params.max:min_max_middle["abs_max"];
          min_max_middle["middle"] = (settings.params.middle !== undefined)?settings.params.middle:data[self._hack_round((len-1)/2)];
        }
        self._prop2settings[prop].values = min_max_middle;

        var values = Object.keys(v2v);
        for(var i = 0, len=values.length; i<len; i++){
          var value = values[i];
          v2v[value] = self._prop2fnc[prop].feature(parseFloat(value));
        }
        v2v[null] = prop2null[prop];
        v2v[undefined] = prop2null[prop];
      }
      self._prop2settings[prop].v2v = v2v;
      return false;
    }
    else{
      console.log(prop, "CACHE");
      return true;
    }
  }

  ChemSpace.prototype._get_min_max = function(data){
    var self = this;
    var i, len;
    var min_max = [];

    var len = data.length;
    data.sort(self._sort_number_ascending);
    min_max.push(data[0]);
    min_max.push(data[data.length-1]);
    return min_max;
  }

  ChemSpace.prototype._get_color_for_value = function(value, min, max, middle, color_scale){
    var self = this;
    if(value === null || value === undefined){
      return self.settings.shapes.circle.fill;
    }

    var colors = [color_scale.min, color_scale.middle, color_scale.max];
    var c1 = colors[0];
    var c2 = colors[colors.length-1];

    if(value > max){
      return c2;
    }

    if(min == max || value < min){
      return c1;
    }

    if(colors.length > 2){
        if(value >= middle){
            min = middle;
            c1 = colors[1];
            c2 = colors[2];
        }
        else{
            max = middle;
            c1 = colors[0];
            c2 = colors[1];
        }
    }
    c1 = self._color_dicts[c1];
    c2 = self._color_dicts[c2];
    var position = (value-min)/(max-min);
    var r = self._hack_round(c1.r+(position*(c2.r-c1.r)));
    var g = self._hack_round(c1.g+(position*(c2.g-c1.g)));
    var b = self._hack_round(c1.b+(position*(c2.b-c1.b)));
    return 'rgb('+r+','+g+','+b+')';
  }

  ChemSpace.prototype._hack_round = function(value){
    var self = this;
      return (0.5 + value) >> 0;
  }

  ChemSpace.prototype._hack_round_float = function(num, dec) { // round number to specified number of decimal places
      var result = Math.round(num*Math.pow(10,dec))/Math.pow(10,dec);
      return result;
  }

  ChemSpace.prototype._sort_number_ascending = function(a, b){
    var self = this;
    return a - b;
  }

  ChemSpace.prototype._sort_number_descending = function(a, b){
    var self = this;
    return b - a;
  }

  ChemSpace.prototype.get_points_by_range = function(range){
    var self = this;
    var key, coords, selected = [];

    for(var i = 0, keys = Object.keys(self.point2coord), len = keys.length; i<len; i++){
      key = keys[i];
      coords = self.point2coord[key];

      if(coords[0] >= range.x[0] && coords[0] <= range.x[1] && coords[1] >= range.y[0] && coords[1] <= range.y[1]){
        // selected = selected.concat(self.data.points[key][self.keys.object_ids]);
        selected.push(key);
      }
    }
    
    return selected;
  }

  ChemSpace.prototype._get_icon_group = function(icon_name, icon_label, options){
    var self = this;
    var icon_group = new Konva.Group({class: "icon"});
    var icon = self.objects_ref.icon.clone({
            data: self.paths_ref[icon_name],
            icon_name: icon_name
    });
    icon.setAttrs(options);

    var icon_overlay = self.objects_ref.icon_overlay.clone({x: options.x, y: options.y, label: icon_label, icon_name: icon_name});
    icon_group.add(icon, icon_overlay);
    return icon_group;
  }

  /************************************
  Compound drawing
  ************************************/
  
  ChemSpace.prototype.draw_compounds = function(point_ids, append){
    var self = this;

    if(append !== undefined && append){
      self.fixed_compounds = self.fixed_compounds.concat(point_ids);
    }
    else{
      self.fixed_compounds = point_ids;
    }
    self._draw_compounds();
  }

  ChemSpace.prototype.erase_compounds = function(){
    var self = this;
    self.fixed_compounds = [];
    self._draw_compounds();
  }

  ChemSpace.prototype._draw_compounds = function(point_ids){
    var self = this;
    
    if(!self.settings.compounds.draw){
        return;
    }

    var key;
    self.compounds_layer.destroyChildren();
    self.compounds_layer.draw();

    if(point_ids === undefined){
      var point_ids = Object.keys(self.point2coord);
      if(point_ids.length > self.settings.compounds.limit){
        if(self.fixed_compounds.length == 0){
          return;
        }
        else{
          point_ids = self.fixed_compounds;
        }
      }
    }

    for(var i = 0, keys=Object.keys(self.point2img), len=keys.length; i<len; i++){
      key = keys[i];

      if(self.point2coord[key] === undefined || point_ids.indexOf(key) === -1){
        self.point2img[key].hide();
      }
      else{
        self.point2img[key].show(); 
      }
    }
    self._draw_compound_structures(point_ids);
  }

  ChemSpace.prototype._prefetch_smiles = async function(ids){
    let self = this;
    let to_get = [];
    let id2object_id = {};
    let object_id2id = {};

    ids.forEach(function(id){
      id2object_id[id] = self.data.points[id][self.keys.object_ids][0];
      object_id2id[self.data.points[id][self.keys.object_ids][0]] = id;
    });

    if(self.data.compounds === undefined){
      self.data.compounds = {};
    }

    ids.forEach(function(id){
      if(self.data.compounds[id] === undefined){
        self.data.compounds[id] = {};
      }
      
      if(self.data.compounds[id][self.keys.smiles] === undefined){
        to_get.push(id2object_id[id]);
      }
    });

    let id2smiles = await self.fetch_smiles(to_get);
    for (const [id, smiles] of Object.entries(id2smiles)) {
      self.data.compounds[object_id2id[id]][self.keys.smiles] = smiles;
    }
    return true;
  }  

  ChemSpace.prototype._draw_compound_structures = async function(point_ids){
    var self = this;
    var structure, key, origin, img;
    var to_get = [];
    var shift = self.settings.compounds.size/2;
    var done = 0;
    let prefetch = await self._prefetch_smiles(point_ids);

    for(var i = 0, keys=point_ids, len=keys.length; i<len; i++){
      key = keys[i];

      if(self.data.compounds[key] !== undefined && self.point2coord[key] !== undefined){
        if(self.point2img[key] === undefined){
          to_get.push(key);
        }
        else{
          self.point2img[key].setAttrs({x: self.point2coord[key][0] - shift, y:self.point2coord[key][1] - shift});
          self.compounds_layer.add(self.point2img[key]);
          self.point2img[key].moveToBottom();
          done++;
        }
      }
    }
    
    function _draw(){
      self._update_visible_compounds();
      self.compounds_layer.draw();
    }

    if(done < to_get.length){
      for(var i = 0, len=to_get.length; i<len; i++){
        var key = to_get[i];
        
        SmilesDrawer.parse(self.data.compounds[key][self.keys.smiles], function (tree) {
          if(self.data.compounds[key].color === undefined){
            self.smilesDrawer.draw(tree, 'chemspacejs-img_cache', 'light', false);
          }
          else{
            var sd_space = $.extend({}, self.settings.compounds.smilesDrawer);
            sd_space.color = self.data.compounds[key].color;
            sd_space.width = self.settings.compounds.size;
            sd_space.height = self.settings.compounds.size;
            current_smilesDrawer = new SmilesDrawer.Drawer(sd_space);
            current_smilesDrawer.draw(tree, 'chemspacejs-img_cache', 'light', false);
          }
          var canvas = document.getElementById("chemspacejs-img_cache");
          var img = new Image();
          img.key = key;
          
          img.onload = function() {
            var image = new Konva.Image({
              image: img,
              width: self.settings.compounds.size,
              height: self.settings.compounds.size,
              listening: false,
              preventDefault: false
            });
            self.point2img[this.key] = image.clone();
            self.point2img[this.key].setAttrs({x: self.point2coord[this.key][0] - shift, y:self.point2coord[this.key][1] - shift});            
            self.compounds_layer.add(self.point2img[this.key]);

            done++;

            if(done == to_get.length){
              _draw();
            }
          };

          img.src = canvas.toDataURL();

        }, function (err) {
            console.log(err);
            done++;

            if(done == to_get.length){
              _draw();
            }
        });
      }
    }
    else{
      _draw();
    }
  }

  ChemSpace.prototype._get_compound_structure = function(point_id){
    var self = this;
    return self.data.compounds[point_id];
  }
  
  ChemSpace.prototype._draw_compound_in_tooltip = async function(point_ids){
    var self = this;
    let prefetch = await self._prefetch_smiles(point_ids);

    for(var i = 0, len=point_ids.length; i<len; i++){
        SmilesDrawer.parse(self.data.compounds[point_ids[i]][self.keys.smiles], function (tree) {
            self.tooltipSmilesDrawer.draw(tree, self.settings.target + "@tooltip_" + i, 'light', false);
        }, function (err) {
        });
    }
  };

  ChemSpace.prototype._draw_target_overlay = function(){
    var self = this;
    var overlay = self.target_element.find(".target_overlay");

    if(overlay.length){
      overlay.fadeIn();
    }
    else{
      overlay = $("<div class='target_overlay'></div>");
      overlay.css({"background-color": "white",
                      "position": "absolute",
                      "top": 0,
                      "left": 0,
                      "right": 0,
                      "bottom": 0,
                      "opacity": 0.5,
                      "width": self.target_element.width(),
                      "height": self.target_element.height(),
                      "z-index": 20
          });
      self.target_element.append(overlay);
    }

    return overlay;
  }

  ChemSpace.prototype._get_average = function(values){
    var sum = 0;
    if(values.length == 1){
      var avg = values[0];
    }
    else{
      for(var i = 0, len=values.length; i<len; i++){
          sum += values[i];
      }
      var avg = sum/values.length;
    }
    return avg;
  };

  ChemSpace.prototype._css_hover = function(elms){
    elms.hover(
      function(){$(this).css({"cursor": "pointer"})},
      function(){$(this).css({})}
    )
  };

  ChemSpace.prototype._percentile_alert = function(values){
    var warning = false;
    if(values.max > 100 || values.max < 0){
      warning = "The maximum percentile must be a number between 0 and 100.";
    }
    else if(values.middle > 100 || values.middle < 0){
      warning = "The middle percentile must be a number between 0 and 100.";
    }
    else if(values.min > 100 || values.min < 0){
      warning = "The minimum percentile must be a number between 0 and 100.";
    }

    if(warning !== false){
      alert(warning);
      return false;
    }
  }

  ChemSpace.prototype._min_max_middle_inputs = function(values_dict, datatype, property){
    var self = this;
    var order = ["min", "middle", "max"];
    var result = $("<div></div>")
      .css({"display": "flex", "justify-content": "space-between", "padding": "8px 0px", "border-bottom": "solid gray 1px", "margin-bottom": 5});

    for(var i = 0, len=order.length; i<len; i++){
      var value = values_dict[order[i]];

      var option = $("<div></div>").css({"width": "30%", "text-align": "center"});
      option.append($("<div>" + self._capitalize(order[i]) + "</div>").css({"font-size": "small", "color": "gray"}));
      var input = $("<div></div>").css({"width": "100%"});
      if(property === "color"){
        input.append($("<input data-type='" + datatype + "' class='option' name='" + order[i] + "' type='color' value='" + self._rgb2hex(value) + "'>")
            .css({"width": "90%", "text-align": "center"})
        );          
      }
      else{
        input.append($("<input data-type='" + datatype + "' class='option' name='" + order[i] + "' type='text' value='" + value + "'>")
            .css({"width": "90%", "text-align": "center"})
        );          
      }
      option.append(input);
      result.append(option);
    }
    return result;
  }

  ChemSpace.prototype._parse_color = function(color){
    var self = this;
    if(color !== undefined){
        if(color.startsWith("#")){
            return self._hex2rgb(color);
        }
        else if(color.startsWith("rgb")){
            return self._rgb2dict(color);
        }
    }
  }

  ChemSpace.prototype._hex2rgb = function(hex) {
    var shorthandRegex = /^#?([a-f\d])([a-f\d])([a-f\d])$/i;
    hex = hex.replace(shorthandRegex, function(m, r, g, b) {
        return r + r + g + g + b + b;
    });

    var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? {
        r: parseInt(result[1], 16),
        g: parseInt(result[2], 16),
        b: parseInt(result[3], 16)
    } : null;
  }

  ChemSpace.prototype._rgb2hex = function(rgb) {
    var self = this;
    var rgb_dict = self._rgb2dict(rgb);
    var hex_color = "#";
    ["r", "g", "b"].forEach(function(key){
        var c = Number(rgb_dict[key]).toString(16);
        if (c.length < 2) {
           c = "0" + c;
        }
        hex_color = hex_color + c;
    });
    return hex_color;
  }

  ChemSpace.prototype._rgb2dict = function(rgb){
    var rgb_split = rgb.split("(")[1].split(")")[0].split(",");
    return {r: parseInt(rgb_split[0]), g: parseInt(rgb_split[1]), b: parseInt(rgb_split[2])};
  }


  ChemSpace.prototype._capitalize = function(str){
    return str.charAt(0).toUpperCase() + str.slice(1);
  }

  ChemSpace.prototype._shuffle = function(arr){
    var newArr = [];
    console.time("SHUFFLE");
    while (arr.length) {

       var randomIndex = Math.floor(Math.random() * arr.length),
           element = arr.splice(randomIndex, 1)

       newArr.push(element[0]);

    }
    console.timeEnd("SHUFFLE");

    return newArr;
  }

  ChemSpace.prototype._translate_objects_to_points = function(objects){
    var self = this;
    var point_ids = [];

    for(var i = 0, len=objects.length; i<len; i++){
      point_ids = point_ids.concat(self.object_id2points[objects[i]]);
    }
    return point_ids;
  }

  ChemSpace.prototype._translate_points_to_objects = function(points){
    var self = this;
    var objects = [];

    for(var i = 0, len=points.length; i<len; i++){
      objects = objects.concat(self.data.points[points[i]][self.keys.object_ids]);
    }
    return objects;
  }

}(jQuery));