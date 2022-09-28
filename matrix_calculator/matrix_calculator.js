matrices = []
matrix_expresison_widgets = []
script_tree = []
method_names = [
    'transpose',
    'determinant',
    'minors',
    'cofactors',
    'adjugate',
    'inverse',
]
function combineMatrices(mat1, mat2, rdiff, cdiff) {
    var br = Math.min(0, rdiff)
    var bc = Math.min(0, cdiff)
    var er = Math.max(mat1.height, mat2.height + rdiff)
    var ec = Math.max(mat1.width, mat2.width + cdiff)
    var matrix_object = createMatrix('', er - br, ec - bc);
    var mat1br = rdiff >= 0 ? 0 : -rdiff;
    var mat1bc = cdiff >= 0 ? 0 : -cdiff;
    var mat2br = rdiff >= 0 ? rdiff : 0;
    var mat2bc = rdiff >= 0 ? cdiff : 0;

    for (var ir = 0; ir < mat1.height; ++ir) {
        for (var ic = 0; ic < mat1.width; ++ic) {
            matrix_object.matrix[mat1br + ir][mat1bc + ic] = mat1.matrix[ir][ic]
        }
    }
    for (var ir = 0; ir < mat2.height; ++ir) {
        for (var ic = 0; ic < mat2.width; ++ic) {
            matrix_object.matrix[mat2br + ir][mat2bc + ic] = mat2.matrix[ir][ic]
        }
    }

    return matrix_object;
}
function getRowFromMatrix(mat, row_index) {
    var matrix_object = createMatrix('', 1, mat.width);
    for (var i = 0; i < mat.width; ++i) {
        matrix_object.matrix[0][i] = mat.matrix[row_index][i];
    }
    return matrix_object;
}
function getColumnFromMatrix(mat, column_index) {
    var matrix_object = createMatrix('', mat.height, 1);
    for (var i = 0; i < mat.height; ++i) {
        matrix_object.matrix[i][0] = mat.matrix[i][column_index];
    }
    return matrix_object;
}
function copyMatrix(dst, src) {
    dst.width = src.width
    dst.height = src.height
    dst.matrix = JSON.parse(JSON.stringify(src.matrix))
}
function cloneMatrix(mat) {
    var matrix_object = createMatrix(mat.name);
    copyMatrix(matrix_object, mat);
    return matrix_object
}
function makeIdentity(size) {
    var matrix_object = createMatrix('I', size, size);
    for (var ir = 0; ir < size; ++ir) {
        for (var ic = 0; ic < size; ++ic) {
            matrix_object.matrix[ir][ic] = ir == ic ? 1 : 0;
        }
    }
    return matrix_object
}
function subMatrix(mat, br, bc, er, ec) {
    var matrix_object = createMatrix(mat.name, er - br, ec - bc);
    for (var ir = br; ir < er; ++ir) {
        for (var ic = bc; ic < ec; ++ic) {
            matrix_object.matrix[ir - br][ic - bc] = mat.matrix[ir][ic];
        }
    }
    return matrix_object;
}
function createMatrix(matrix_name, height = 2, width = 2) {
    var matrix = new Array(height);
    for (var i = 0; i < matrix.length; ++i) {
        matrix[i] = new Array(width).fill(0);
    }
    return {
        name: matrix_name,
        type: 'matrix',
        height: height,
        width: width,
        matrix: matrix,
    }; // row-major
}
function addMatrices(a, b) {
    var r = a.height > b.height ? a.height : b.height;
    var c = a.width > b.width ? a.width : b.width;
    var matrix_object = createMatrix('?', r, c)
    reserveMatrix(a, r, c);
    reserveMatrix(b, r, c);
    for (var ir = 0; ir < r; ++ir) {
        for (var ic = 0; ic < c; ++ic) {
            matrix_object.matrix[ir][ic] = a.matrix[ir][ic] + b.matrix[ir][ic]
        }
    }
    return matrix_object
}
function subMatrices(a, b) {
    var r = a.height > b.height ? a.height : b.height;
    var c = a.width > b.width ? a.width : b.width;
    var matrix_object = createMatrix('?', r, c)
    reserveMatrix(a, r, c);
    reserveMatrix(b, r, c);
    for (var ir = 0; ir < r; ++ir) {
        for (var ic = 0; ic < c; ++ic) {
            matrix_object.matrix[ir][ic] = a.matrix[ir][ic] - b.matrix[ir][ic]
        }
    }
    return matrix_object
}
function multMatrices(a, b) {
    if (a == null && b == null) throw 'multMatrices() : a & b cannot be null'
    if (a == null) throw 'multMatrices() : a cannot be null'
    if (b == null) throw 'multMatrices() : b cannot be null'
    if (a.width == b.height) {
        var matrix_object = createMatrix('', a.height, b.width);
        for (var ir = 0; ir < a.height; ++ir) {
            for (var ic = 0; ic < b.width; ++ic) {
                var sum = 0;
                for (var j = 0; j < a.width; ++j) {
                    sum += a.matrix[ir][j] * b.matrix[j][ic]
                }
                matrix_object.matrix[ir][ic] = sum
            }
        }
        return matrix_object
    }
    else {
        throw `cannot multiply matrices of size (${a.height}×${a.width}),(${b.height}×${b.width}).`
        // return null
    }
}
function scaleMatrix(k, mat) {
    var matrix_object = createMatrix('', mat.height, mat.width)
    for (var ir = 0; ir < mat.height; ++ir) {
        for (var ic = 0; ic < mat.width; ++ic) {
            matrix_object.matrix[ir][ic] = k * mat.matrix[ir][ic]
        }
    }
    return matrix_object
}
function makeMinor(mat, r, c) {
    var matrix_object = createMatrix('', mat.height - 1, mat.width - 1);
    var idstr = 0;
    for (var ir = 0; ir < mat.height; ++ir) {
        if (ir == r) {
            continue;
        }
        var idstc = 0;
        for (var ic = 0; ic < mat.width; ++ic) {
            if (c == ic)
                continue
            matrix_object.matrix[idstr][idstc++] = mat.matrix[ir][ic];
        }
        ++idstr;
    }
    return matrix_object;
}
function makeMinors(mat) {
    if (mat.width != mat.height)
        throw 'makeMinors() : non-square matrices do not have minors'
    var matrix_object = createMatrix('', mat.height, mat.width);
    for (var ir = 0; ir < mat.height; ++ir) {
        for (var ic = 0; ic < mat.width; ++ic) {
            var minor = makeMinor(mat, ir, ic)
            matrix_object.matrix[ir][ic] = calcDeterminant(minor)
        }
    }
    return matrix_object;
}
function makeCofactors(mat) {
    if (mat.width != mat.height)
        throw 'makeCofactors() : non-square matrices do not have cofactors'
    var matrix_object = createMatrix('', mat.height, mat.width);
    for (var ir = 0; ir < mat.height; ++ir) {
        var sign = ir % 2 == 0 ? 1 : -1;
        for (var ic = 0; ic < mat.width; ++ic) {
            var minor = makeMinor(mat, ir, ic)
            matrix_object.matrix[ir][ic] = sign * calcDeterminant(minor)
            sign *= -1;
        }
    }
    return matrix_object;
}
function makeInverse(mat) {
    if (mat.width != mat.height)
        throw 'makeInverse() : non-square matrices do not have inverses'
    var determinant = calcDeterminant(mat);
    var epsilon = 0.000000001;
    if (-epsilon < determinant && determinant < epsilon)
        throw 'makeInverse() : singular matrices do not have inverses (determinant is zero)'
    var adjugate = makeTranspose(makeCofactors(mat))
    return scaleMatrix(1 / determinant, adjugate);
}
function calcDeterminant(mat) {
    if (mat.width != mat.height) {
        throw `calcDeterminant() : non-square matrices do not have determinants. (${mat.height}×${mat.width})`
    }
    if (mat.width == 1 && mat.height == 1) {
        return mat.matrix[0][0]
    }
    var sign = 1;
    var determinant = 0;
    for (var ic = 0; ic < mat.width; ++ic) {
        var a = mat.matrix[0][ic]
        var minor = makeMinor(mat, 0, ic)
        var minor_determinant = calcDeterminant(minor)
        var c = sign * a * minor_determinant
        determinant += c
        sign *= -1
    }
    return determinant
}
function makeTranspose(mat) {
    var matrix_object = createMatrix('', mat.width, mat.height)
    for (var ir = 0; ir < mat.height; ++ir) {
        for (var ic = 0; ic < mat.width; ++ic) {
            matrix_object.matrix[ic][ir] = mat.matrix[ir][ic]
        }
    }
    return matrix_object
}
function findMatrixObjectByName(matrix_name) {
    return matrices.find(element => element.name == matrix_name)
}
function reserveMatrix(matrix_object, rows_count, cols_count) {
    if (getMatrixRowsReserve(matrix_object) < rows_count) {
        var lacks = rows_count - getMatrixRowsReserve(matrix_object)
        var tail = Array(getMatrixColsReserve(matrix_object)).fill(0)
        matrix_object.matrix = matrix_object.matrix.concat(Array(lacks).fill(tail))
    }
    if (getMatrixColsReserve(matrix_object) < cols_count) {
        var lacks = cols_count - getMatrixColsReserve(matrix_object)
        for (var ir = 0; ir < matrix_object.matrix.length; ++ir) {
            matrix_object.matrix[ir] = matrix_object.matrix[ir].concat(Array(lacks).fill(0));
        }
    }
}
function resizeMatrix(matrix_object, rows_count, cols_count) {
    reserveMatrix(matrix_object, rows_count, cols_count);
    matrix_object.height = rows_count;
    matrix_object.width = cols_count;
}
function getMatrixColsReserve(matrix_object) {
    return matrix_object.matrix[0].length
}
function getMatrixRowsReserve(matrix_object) {
    return matrix_object.matrix.length
}
var fpCreateTable = (r, c, fpTdContentCreator) => {
    var table = document.createElement('table');
    var tbody = document.createElement('tbody');
    table.appendChild(tbody)
    for (var ir = 0; ir < r; ++ir) {
        var tr = document.createElement('tr');
        for (var ic = 0; ic < c; ++ic) {
            var td = document.createElement('td');
            var element_td_content = fpTdContentCreator(ir, ic, ir * c + ic);
            td.appendChild(element_td_content);
            tr.appendChild(td);
        }
        tbody.appendChild(tr);
    }
    return table;
}
function isWhitespaceChar(c) {
    if (c == undefined) return false
    return /\s/.test(c)
}
function isOperatorChar(c) {
    if (c == undefined) return false
    if (c == '=' || c == '+' || c == '-' || c == '*' || c == '^' || c == '/' || c == '@' || c == ';' || c == '(' || c == ')' || c == '[' || c == ']' || c == '{' || c == '}' || c == ',') return true
    return false
}
function isNameChar(c) {
    if (c == undefined) return false
    return !isWhitespaceChar(c) && !isOperatorChar(c)
}
function markWhitespace(str, from) {
    var i = from
    while (isWhitespaceChar(str[i++]));
    return [from, i - 1, str.substr(from, i - from - 1), 'whitespace']
}
function markName(str, from) {
    var i = from
    while (isNameChar(str[i++]));
    return [from, i - 1, str.substr(from, i - from - 1), 'name']
}
function markMethod(str, from) {
    for (var i = 0; i < method_names.length; ++i) {
        var method_name = method_names[i]
        var substr = str.substr(from, method_name.length)
        if (substr == method_name) {
            return [from, from + method_name.length, method_name, 'method'];
        }
    }
    return [from, from, null, 'method'];
}
function markOperator(str, from) {
    if (isOperatorChar(str[from]))
        return [from, from + 1, str[from], 'operator']
    else
        return [from, from, null, 'operator']
}
function markNumber(str, from) {
    var float_regex = /[\+\-]?((\d+(\.\d*)?)|(\.\d+))([eE][\+\-]?\d*)?/g
    float_regex.lastIndex = from
    var ary = float_regex.exec(str);
    if (ary == null || ary.index != from) {
        return [from, from, NaN, 'number']
    }
    var number = parseFloat(ary[0])
    if (isNaN(number)) {
        console.warn('markNumber() : parseFloat() failed while having found a regex match')
        return [from, from, NaN, 'number']
    }
    return [from, from + ary[0].length, number, 'number']
}
function markSomething(str, from) {
    var mark;
    mark = markWhitespace(str, from)
    if (mark[0] != mark[1])
        return mark
    mark = markNumber(str, from)
    if (mark[0] != mark[1])
        return mark
    mark = markMethod(str, from)
    if (mark[0] != mark[1])
        return mark
    mark = markOperator(str, from)
    if (mark[0] != mark[1])
        return mark
    mark = markName(str, from)
    if (mark[0] != mark[1])
        return mark
    return null
}
function markAll(str) {
    if (str == null)
        return null;
    i = 0
    marks = []
    while (i != str.length) {
        var mark = markSomething(str, i)
        if (mark[0] == mark[1])
            return null; // fail
        marks.push(mark)
        i = mark[1]
    }
    return marks
}
function onLoad(event) {
    var template_matrix_input = document.getElementById("template_matrix_input")
    var template_matrix_display_entry = document.getElementById("template_matrix_display_entry")
    var template_number_display_entry = document.getElementById("template_number_display_entry")
    var element_matrix_area = document.getElementById("matrix_area")
    var element_add_widget_button = document.getElementById("add_widget_button")
    var element_writer_section_textarea = document.getElementsByClassName('writer_section_textarea')[0]
    var element_matrix_display = document.getElementById('matrix_display')
    var element_menu_exporter = document.getElementById('menu_exporter')
    element_writer_section_textarea.addEventListener('keyup', (event) => {
        if (event.target == element_writer_section_textarea) {
            recalcMatrix()
        }
    })
    loadFromState(loadStateFromLocalStorage());
    function evalStatement(marks) {
        if (marks == null) throw 'evalStatement() : marks cannot be null'
        nodes = marks.filter(mark => mark[3] != 'whitespace').map(mark => {
            return {
                l: null,
                r: null,
                mark: mark,
                val: null,
            }
        })
        if (nodes.length > 0) {
            root = treeficate(nodes)
            return evalNode(root)
        }
        return {
            val: null,
            val_type: 'empty',
            node: null,
        }
    }
    function renewMatrixList(marks) {
        var marks = marks.filter(mark => mark[3] == 'name')
        for (var i = 0; i < marks.length; ++i) {
            var mark = marks[i];
            var name = mark[2];
            var mo = matrices.findIndex(m => m.name == name)
            if (mo == -1) {
                matrices.push(createMatrix(name))
            }
        }
        for (var i = matrices.length - 1; i >= 0; --i) {
            var mo = matrices[i]
            if (marks.findIndex(mark => mark.name == mo.name) == -1 && !mo.input_node) {
                matrices.splice(i, 1);
            }
        }
    }
    function evalScript() {
        str = element_writer_section_textarea.value // something like '   C = A + B   ...'
        if (str == '')
            return;
        var marks = markAll(str)
        renewMatrixList(marks)
        var statements = divideStatements(marks)
        var results = statements.map(statement => evalStatement(statement));
        element_matrix_display.innerText = ''
        displayResults(statements, results);
    }
    function recalcMatrix() {
        evalScript()
        var updates = document.getElementsByClassName('update_table_on_recalc')
        for (var i = 0; i < updates.length; ++i) {
            var name = updates[i].getAttribute("matrix_name")
            var mat = findMatrixObjectByName(name)
            if (mat == null)
                continue;
            var table = fpCreateTable(mat.height, mat.width, (ir, ic, i) => {
                var label = document.createElement('span');
                var element = mat.matrix[ir][ic];
                label.innerHTML = element;
                return label;
            })
            updates[i].innerHTML = ""
            updates[i].appendChild(table)
        }
        saveStateToLocalStorage()
    }
    function createMatrixInputNode(ir, ic, i, matrix_object) {
        var input = document.createElement('input');
        input.setAttribute('type', 'number');
        input.addEventListener('input', () => {
            matrix_object.matrix[ir][ic] = input.valueAsNumber
            recalcMatrix()
        });
        input.classList.add('matrix_element_input')
        input.value = matrix_object.matrix[ir][ic]
        return input
    }
    function establishMatrixInput(matrix_input_node, matrix_object) {
        var element_input_size_rows = matrix_input_node.getElementsByClassName("matrix_input_size_rows")[0]
        var element_input_size_cols = matrix_input_node.getElementsByClassName("matrix_input_size_cols")[0]
        var element_input_matrix_input_elements = matrix_input_node.getElementsByClassName("matrix_input_elements")[0]
        var rows_count = element_input_size_rows.valueAsNumber
        var cols_count = element_input_size_cols.valueAsNumber
        resizeMatrix(matrix_object, rows_count, cols_count)
        var table = fpCreateTable(rows_count, cols_count, (ir, ic, i) => {
            var node = createMatrixInputNode(ir, ic, i, matrix_object)
            return node;
        })
        table.classList.add('matrix_input_element_table')
        element_input_matrix_input_elements.innerHTML = ""
        element_input_matrix_input_elements.appendChild(table)
        matrix_object.input_node = matrix_input_node
    }
    function createMatrixInputWidget(matrix_object) {
        var div = document.createElement('div')
        div.appendChild(template_matrix_input.content.cloneNode(true))
        var widget = div.children[0]
        widget.getElementsByClassName('matrix_input_remove_button')[0].addEventListener('click', () => {
            matrices.splice(matrices.indexOf(matrix_object), 1);
            widget.parentElement.removeChild(widget)
        })
        widget.getElementsByClassName('matrix_input_name')[0].textContent = matrix_object.name
        var element_input_size_rows = widget.getElementsByClassName("matrix_input_size_rows")[0]
        var element_input_size_cols = widget.getElementsByClassName("matrix_input_size_cols")[0]
        element_input_size_rows.value = matrix_object.height
        element_input_size_rows.addEventListener('input', () => {
            matrix_object.height = element_input_size_rows.valueAsNumber
            establishMatrixInput(widget, matrix_object)
            recalcMatrix()
        })
        element_input_size_cols.value = matrix_object.width;
        element_input_size_cols.addEventListener('input', () => {
            matrix_object.width = element_input_size_cols.valueAsNumber
            establishMatrixInput(widget, matrix_object)
            recalcMatrix()
        })
        establishMatrixInput(widget, matrix_object)
        return widget
    }
    recalcMatrix()
    window.onclick = function (event) {
        if (event.target.matches('.export_button > svg')) {
            event.target.parentElement.appendChild(element_menu_exporter)
            element_menu_exporter.classList.add('show')
        }
        else {
            if (element_menu_exporter.classList.contains('show')) {
                element_menu_exporter.classList.remove('show');
            }
        }
        if (event.target.matches('#a_json_row_major')) {
            var widget = event.target.parentElement.parentElement.parentElement;
            var element_matrix_output = widget.getElementsByClassName('matrix_output')[0]
            var name = element_matrix_output.getAttribute("matrix_name")
            var matrix = matrices.find(m => m.name == name)
            if (matrix == null)
                matrix = widget.targetMatrix
            var to_stringify = subMatrix(matrix, 0, 0, matrix.height, matrix.width);
            navigator.clipboard.writeText(JSON.stringify(to_stringify.matrix))
            event.preventDefault()
        }
        if (event.target.matches('#a_json_col_major')) {
            var widget = element_menu_exporter.parentElement.parentElement;
            var element_matrix_output = widget.getElementsByClassName('matrix_output')[0]
            var name = element_matrix_output.getAttribute("matrix_name")
            var matrix = matrices.find(m => m.name == name)
            if (matrix == null)
                matrix = widget.targetMatrix
            navigator.clipboard.writeText(JSON.stringify(makeTranspose(subMatrix(matrix, 0, 0, matrix.height, matrix.width)).matrix))
            event.preventDefault()
        }
        if (event.target.matches('#a_katex')) {
            var widget = element_menu_exporter.parentElement.parentElement;
            var element_matrix_output = widget.getElementsByClassName('matrix_output')[0]
            var name = element_matrix_output.getAttribute("matrix_name")
            var matrix = matrices.find(m => m.name == name)
            if (matrix == null)
                matrix = widget.targetMatrix
            var str = '\\begin{bmatrix}'
            for (var ir = 0; ir < matrix.height; ++ir) {
                for (var ic = 0; ic < matrix.width; ++ic) {
                    var element = matrix.matrix[ir][ic]
                    str += element.toString()
                    if (ic < matrix.width - 1)
                        str += '&'
                    else
                        str += '\\\\'
                }
            }
            str += '\\end{bmatrix}'
            navigator.clipboard.writeText(str)
            event.preventDefault()
        }
        if (event.target.matches('#a_c_row_major')) {
            var widget = element_menu_exporter.parentElement.parentElement;
            var element_matrix_output = widget.getElementsByClassName('matrix_output')[0]
            var name = element_matrix_output.getAttribute("matrix_name")
            var matrix = matrices.find(m => m.name == name)
            if (matrix == null)
                matrix = widget.targetMatrix
            var str = JSON.stringify(matrix.matrix).replaceAll('[', '{').replaceAll(']', '}')
            navigator.clipboard.writeText(str)
            event.preventDefault()
        }
        if (event.target.matches('#a_c_col_major')) {
            var widget = element_menu_exporter.parentElement.parentElement;
            var element_matrix_output = widget.getElementsByClassName('matrix_output')[0]
            var name = element_matrix_output.getAttribute("matrix_name")
            var matrix = makeTranspose(matrices.find(m => m.name == name))
            if (matrix == null)
                matrix = widget.targetMatrix
            var str = JSON.stringify(matrix.matrix).replaceAll('[', '{').replaceAll(']', '}')
            navigator.clipboard.writeText(str)
            event.preventDefault()
        }
        if (event.target.matches('#add_widget_button')) {
            var names = matrices.map(matrix_object => matrix_object.name)

            function nextChar(char) {
                if (char == 'H') return 'J';
                if (char == 'N') return 'P';
                return String.fromCharCode(char.charCodeAt(0) + 1);
            }
            function nextName(name) {
                var new_name_ary = []
                var carry = 1;
                for (var i = name.length - 1; i >= 0; --i) {
                    if (carry == 0) {
                        new_name_ary.push(name[i])
                    }
                    else {
                        if (name[i] == 'Z') {
                            carry = 1
                            new_name_ary.push('A')
                        }
                        else {
                            carry = 0
                            new_name_ary.push(nextChar(name[i]));
                        }
                    }
                }
                if (carry == 1) {
                    new_name_ary.push('A')
                }
                return new_name_ary.reverse().join('')
            }
            var new_name = 'A';
            while (names.indexOf(new_name) != -1) {
                new_name = nextName(new_name);
            }
            var matrix_object = createMatrix(new_name)
            matrices.push(matrix_object)
            element_matrix_area.insertBefore(createMatrixInputWidget(matrix_object), element_add_widget_button)
            event.preventDefault()
        }
        if (event.target.matches('.top_button_example')) {
            for (var i = 0; i < matrices.length; ++i) {
                var input_node = matrices[i].input_node;
                if( input_node )
                {
                    input_node.parentElement.removeChild(input_node);
                }
            }
            matrices.splice(0, matrices.length);

            var sstate = {
                "matrices": "[{\"name\":\"A\",\"type\":\"matrix\",\"width\":3,\"height\":3,\"matrix\":[[3,0,2,0],[2,5,5,3],[0,1,1,3],[5,-18,1,0]],\"input_node\":{}},{\"name\":\"B\",\"type\":\"matrix\",\"height\":3,\"width\":3,\"matrix\":[[1,2,3],[4,5,6],[7,8,9]],\"input_node\":{}},{\"name\":\"C\",\"type\":\"matrix\",\"height\":2,\"width\":2,\"matrix\":[[1,2,0],[3,4,0],[0,0,0]],\"input_node\":{}},{\"name\":\"D\",\"type\":\"matrix\",\"height\":3,\"width\":3,\"matrix\":[[1,1,1],[1,1,1],[1,1,2]],\"input_node\":{}}]",
                "expressions": "A{0};\nA(0);\ntranspose(A);\ntranspose({1,2,3});\n\ninverse(A);\nminors(A);\ncofactors(A);\ninverse(A);\n\nA(0,0);\nA(2,0);\nA(1,2);\n(1,2,3);\n{1,2,3};\n{(1,2,3),(4,5,6),(7,8,9)};\n({1,2,3},{4,5,6},{7,8,9});\n({(1,2,3),(4,5,6),(7,8,9)},({1,2,3},{4,5,6},{7,8,9}));\n\nC@D;\nD@C;\n"
            }
            loadFromState(sstate)
            recalcMatrix();
        }
    }
    function divideStatements(marks) {
        var statements = []
        var statement = []
        for (var i = 0; i < marks.length; ++i) {
            var mark = marks[i]
            if (mark[2] == ';') {
                statements.push(statement)
                statement = []
                continue;
            }
            else {
                statement.push(marks[i])
            }
        }
        return statements
    }
    function treeficate(nodes) {
        var to_fold = []
        var depth = 0;
        for (var i = 0; i < nodes.length; ++i) {
            var node = nodes[i]
            // console.log(node.mark, depth);
            if (node.mark[2] == '(' || node.mark[2] == '[' || node.mark[2] == '{') {
                ++depth;
                to_fold.push({
                    node: node,
                    depth: depth
                })
            }
            else if (node.mark[2] == ')' || node.mark[2] == ']' || node.mark[2] == '}') {
                to_fold.push({
                    node: node,
                    depth: depth
                })
                --depth;
            }
            else if (node.mark[3] == 'method') {
                to_fold.push({
                    node: node,
                    depth: depth
                })
            }
            else if (node.mark[3] == 'operator') {
                to_fold.push({
                    node: node,
                    depth: depth
                })
            }
            else if (node.mark[3] == 'name') {
                if (i < nodes.length - 1 && nodes[i + 1].mark[3] == 'number') {
                    to_fold.push({
                        node: node,
                        depth: depth
                    })
                }
            }
        }
        var optier =
        {
            '=': 0,
            '(': 1,
            '[': 1,
            '{': 1,
            ')': 2,
            ']': 2,
            '}': 2,
            ',': 3,
            '+': 4,
            '-': 4,
            '*': 5,
            '/': 5,
            '@': 5,
            '^': 6,
            'method': 7,
        };
        fpHash = function (a) {
            var str = a.node.mark[2]
            var type = a.node.mark[3]
            if (type == 'operator') {
                return a.depth * 1000 + optier[str]
            }
            else if (type == 'method') {
                return a.depth * 1000 + optier['method']
            }
            throw 'type must be operator|method'
        }
        to_fold.sort((a, b) => {
            return fpHash(a) - fpHash(b)
        })
        // console.log(to_fold);
        var close_bracket_queue = [];
        for (var d = to_fold.length - 1; d >= 0; --d) {
            var fold = to_fold[d]
            var opnode = fold.node
            var index = nodes.indexOf(opnode)
            if (opnode.mark[3] == 'operator') {
                if (opnode.mark[2] == '=' || opnode.mark[2] == '+' || opnode.mark[2] == '-' || opnode.mark[2] == '*' || opnode.mark[2] == '/' || opnode.mark[2] == '@' || opnode.mark[2] == '^' || opnode.mark[2] == ',') {
                    opnode.r = nodes.splice(index + 1, 1)[0];
                    opnode.l = nodes.splice(index - 1, 1)[0];
                }
                else if (opnode.mark[2] == ']') {
                    close_bracket_queue.push(opnode);
                    nodes.splice(index, 1);
                    opnode.l = nodes.splice(index - 1, 1)[0];
                }
                else if (opnode.mark[2] == '[') {
                    opnode.l = nodes.splice(index - 1, 1)[0];
                    opnode.r = close_bracket_queue.splice(0, 1)[0]
                }
                else if (opnode.mark[2] == ')') {
                    close_bracket_queue.push(opnode);
                    nodes.splice(index, 1);
                    opnode.l = nodes.splice(index - 1, 1)[0];
                }
                else if (opnode.mark[2] == '(') {
                    if (index > 0 && nodes[index - 1].mark[3] == 'name') {
                        opnode.l = nodes.splice(index - 1, 1)[0];
                        opnode.r = close_bracket_queue.splice(0, 1)[0]
                    }
                    else {
                        opnode.r = close_bracket_queue.splice(0, 1)[0]
                    }
                }
                else if (opnode.mark[2] == '}') {
                    close_bracket_queue.push(opnode);
                    nodes.splice(index, 1);
                    opnode.l = nodes.splice(index - 1, 1)[0];
                }
                else if (opnode.mark[2] == '{') {
                    if (index > 0 && nodes[index - 1].mark[3] == 'name') {
                        opnode.l = nodes.splice(index - 1, 1)[0];
                        opnode.r = close_bracket_queue.splice(0, 1)[0]
                    }
                    else {
                        opnode.r = close_bracket_queue.splice(0, 1)[0]
                    }
                }
            }
            else if (opnode.mark[3] == 'method') {
                opnode.r = nodes.splice(index + 1, 1)[0];
            }
            else if (opnode.mark[3] == 'name') {
                opnode.r = nodes.splice(index + 1, 1)[0];
            }
        }
        return nodes[0]
    }
    function evalNode(node) {
        if (node == null)
            throw 'evalNode() : node is null'
        if (node.mark == null)
            throw 'evalNode() : node.mark is null'

        if (node.mark[2] == '=') {
            var le = evalNode(node.l)
            var re = evalNode(node.r)
            if (le.val_type == 'matrix' && re.val_type == 'matrix') {
                if (le.val.input_node != null)
                    throw 'evalNode() : cannot assign to predefined matrix'
                copyMatrix(le.val, re.val)
                return {
                    val_type: 'matrix',
                    val: le.val,
                    node: node
                }
            }
            else if (le.val_type == 'matrix' && re.val_type == 'column') {
                if (le.val.input_node != null)
                    throw 'evalNode() : cannot assign to predefined matrix'
                le.val.matrix = re.val.map(e => [e])
                le.val.height = re.val.length
                le.val.width = 1
                return {
                    val_type: 'matrix',
                    val: le.val,
                    node: node
                }
            }
            else if (le.val_type == 'matrix' && re.val_type == 'columns') {
                if (le.val.input_node != null)
                    throw 'evalNode() : cannot assign to predefined matrix'
                var height = re.val[0].length;
                var width = re.val.length;
                resizeMatrix(le.val, height, width);
                for (var ir = 0; ir < height; ++ir) {
                    for (var ic = 0; ic < width; ++ic) {
                        le.val.matrix[ir][ic] = re.val[ic][ir];
                    }
                }
                return {
                    val_type: 'matrix',
                    val: le.val,
                    node: node
                }
            }
            else if (le.val_type == 'matrix' && re.val_type == 'row') {
                if (le.val.input_node != null)
                    throw 'evalNode() : cannot assign to predefined matrix'
                le.val.matrix = [re.val];
                le.val.height = 1
                le.val.width = re.val.length;
                return {
                    val_type: 'matrix',
                    val: le.val,
                    node: node
                }
            }
            else if (le.val_type == 'matrix' && re.val_type == 'rows') {
                if (le.val.input_node != null)
                    throw 'evalNode() : cannot assign to predefined matrix'
                le.val.matrix = re.val
                le.val.height = re.val.length
                le.val.width = re.val[0].length
                return {
                    val_type: 'matrix',
                    val: le.val,
                    node: node
                }
            }
        }
        else if (node.mark[2] == '+') {
            var le = evalNode(node.l)
            var re = evalNode(node.r)
            if (le.val_type == 'matrix' && re.val_type == 'matrix') {
                return {
                    val_type: 'matrix',
                    val: addMatrices(le.val, re.val),
                    node: node
                }
            }
            else if (le.val_type == 'number' && re.val_type == 'number') {
                return {
                    val_type: 'number',
                    val: le.val + re.val,
                    node: node
                }
            }
        }
        else if (node.mark[2] == '-') {
            var le = evalNode(node.l)
            var re = evalNode(node.r)
            if (le.val_type == 'matrix' && re.val_type == 'matrix') {
                return {
                    val_type: 'matrix',
                    val: subMatrices(le.val, re.val),
                    node: node
                }
            }
            else if (le.val_type == 'number' && re.val_type == 'number') {
                return {
                    val_type: 'number',
                    val: le.val - re.val,
                    node: node
                }
            }
        }
        else if (node.mark[2] == '*') {
            var le = evalNode(node.l)
            var re = evalNode(node.r)
            if (le.val_type == 'matrix' && re.val_type == 'matrix') {
                return {
                    val_type: 'matrix',
                    val: multMatrices(le.val, re.val),
                    node: node
                }
            }
            else if (le.val_type == 'number' && re.val_type == 'matrix') {
                return {
                    val_type: 'matrix',
                    val: scaleMatrix(le.val, re.val),
                    node: node
                }
            }
            else if (le.val_type == 'matrix' && re.val_type == 'number') {
                return {
                    val_type: 'matrix',
                    val: scaleMatrix(re.val, le.val),
                    node: node
                }
            }
            else if (le.val_type == 'number' && re.val_type == 'number') {
                return {
                    val_type: 'number',
                    val: le.val * re.val,
                    node: node
                }
            }
        }
        else if (node.mark[2] == '/') {
            var le = evalNode(node.l)
            var re = evalNode(node.r)
            if (le.val_type == 'matrix' && re.val_type == 'matrix') {
                return {
                    val_type: 'matrix',
                    val: multMatrices(le.val, makeInverse(re.val)),
                    node: node
                }
            }
            else if (le.val_type == 'number' && re.val_type == 'matrix') {
                // throw 'evalNode() : syntax_error. dividing number by matrix is not defined.'
                return {
                    val_type: 'matrix',
                    val: scaleMatrix(le.val, makeInverse(re.val)),
                    node: node
                }
            }
            else if (le.val_type == 'matrix' && re.val_type == 'number') {
                return {
                    val_type: 'matrix',
                    val: scaleMatrix(1 / re.val, le.val),
                    node: node
                }
            }
            else if (le.val_type == 'number' && re.val_type == 'number') {
                return {
                    val_type: 'number',
                    val: le.val / re.val,
                    node: node
                }
            }
        }
        else if (node.mark[2] == '@') {
            var le = evalNode(node.l)
            var re = evalNode(node.r)
            if (le.val_type == 'matrix' && re.val_type == 'matrix') {
                var lhs = le.val;
                var rhs = re.val;
                var matrix_object = createMatrix('');
                for (var ir = 0; ir < lhs.height; ++ir) {
                    for (var ic = 0; ic < lhs.width; ++ic) {
                        var scaled = scaleMatrix(lhs.matrix[ir][ic], rhs)
                        matrix_object = combineMatrices(matrix_object, scaled, ir * rhs.height, ic * rhs.width);
                    }
                }
                return {
                    val_type: 'matrix',
                    val: matrix_object,
                    node: node
                }
            }
        }
        else if (node.mark[2] == '^') {
            var le = evalNode(node.l)
            var re = evalNode(node.r)
            if (le.val_type == 'matrix' && re.val_type == 'matrix') {
                throw 'evalNode() : syntax_error. cannot use matrix as exponent.'
            }
            else if (le.val_type == 'number' && re.val_type == 'matrix') {
                throw 'evalNode() : syntax_error. cannot use matrix as exponent.'
            }
            else if (le.val_type == 'matrix' && re.val_type == 'number') {
                if (le.val.width != le.val.height)
                    throw `evalNode() : syntax_error. cannot take powers of non-square matrices. size is (${le.val.height}x${le.val.width})`
                var m = makeIdentity(le.val.height);
                for (var i = 0; i < re.val; ++i) {
                    m = multMatrices(m, le.val);
                }
                return {
                    val_type: 'matrix',
                    val: m,
                    node: node
                }
            }
            else if (le.val_type == 'number' && re.val_type == 'number') {
                return {
                    val_type: 'number',
                    val: Math.pow(le.val, re.val),
                    node: node
                }
            }
        }
        else if (node.mark[2] == '(') {
            if (node.l != null && node.l.mark[3] == 'name') {
                var le = evalNode(node.l);
                var re = evalNode(node.r);
                if (le.val_type == 'matrix', re.val_type == 'number') {
                    return {
                        val_type: 'matrix',
                        val: getColumnFromMatrix(le.val, re.val),
                        node: node
                    }
                }
                else if (le.val_type == 'matrix', re.val_type == 'numbers') {
                    return {
                        val_type: 'number',
                        val: le.val.matrix[re.val[0]][re.val[1]],
                        node: node
                    }
                }
            }
            else if (node.l == null) {
                var re = evalNode(node.r);
                if (re.val_type == 'number') {
                    return {
                        val_type: 'number',
                        val: re.val,
                        node: node
                    }
                }
                else if (re.val_type == 'numbers') {
                    var mat = createMatrix('', re.val.length, 1);
                    for (var i = 0; i < re.val.length; ++i) {
                        mat.matrix[i][0] = re.val[i];
                    }
                    return {
                        val_type: 'matrix',
                        val: mat,
                        node: node
                    }
                }
                else if (re.val_type == 'matrix') {
                    return {
                        val_type: 'matrix',
                        val: re.val,
                        node: node
                    }
                }
                else if (re.val_type == 'matrices') {
                    var height = 0;
                    var width = 0;
                    for (var i = 0; i < re.val.length; ++i) {
                        var mat = re.val[i];
                        width = mat.width > width ? mat.width : width;
                        height = height + mat.height;
                    }
                    var matrix_object = createMatrix('', height, width);
                    var br = 0;
                    for (var i = 0; i < re.val.length; ++i) {
                        var mat = re.val[i];
                        var er = br + mat.height;
                        var ec = matrix_object.width;
                        for (var ir = br; ir < er; ++ir) {
                            for (var ic = 0; ic < ec; ++ic) {
                                matrix_object.matrix[ir][ic] = mat.matrix[ir - br][ic];
                            }
                        }
                        br = er;
                    }
                    return {
                        val_type: 'matrix',
                        val: matrix_object,
                        node: node
                    }
                }
            }
            else if (re.val_type == 'matrix') {
                return {
                    val_type: 'matrix',
                    val: re.val,
                    node: node
                }
            }
            else if (re.val_type == 'matrices') {
                var height = 0;
                var width = 0;
                for (var i = 0; i < re.val.length; ++i) {
                    var mat = re.val[i];
                    height = mat.height > height ? mat.height : height;
                    width = width + mat.width;
                }
                var matrix_object = createMatrix('', height, width);
                var bc = 0;
                for (var i = 0; i < re.val.length; ++i) {
                    var mat = re.val[i];
                    var er = mat.height;
                    var ec = bc + mat.width;
                    for (var ir = 0; ir < er; ++ir) {
                        for (var ic = bc; ic < ec; ++ic) {
                            matrix_object.matrix[ir][ic] = mat.matrix[ir][ic - bc];
                        }
                    }
                    bc = ec;
                }
                return {
                    val_type: 'matrix',
                    val: matrix_object,
                    node: node
                }
            }
        }
        else if (node.mark[2] == ')') {
            return evalNode(node.l);
        }
        else if (node.mark[2] == '{') {
            if (node.l) {
                if (node.l.mark[3] == 'name') {
                    var le = evalNode(node.l);
                    var re = evalNode(node.r);
                    if (le.val_type == 'matrix', re.val_type == 'number') {
                        return {
                            val_type: 'matrix',
                            val: getRowFromMatrix(le.val, re.val),
                            node: node
                        }
                    }
                }
            } else {
                var re = evalNode(node.r);
                if (re.val_type == 'number') {
                    return {
                        val_type: 'number',
                        val: re.val,
                        node: node
                    }
                }
                else if (re.val_type == 'numbers') {
                    var mat = createMatrix('', 1, re.val.length);
                    for (var i = 0; i < re.val.length; ++i) {
                        mat.matrix[0][i] = re.val[i];
                    }
                    return {
                        val_type: 'matrix',
                        val: mat,
                        node: node
                    }
                }
                else if (re.val_type == 'matrix') {
                    return {
                        val_type: 'matrix',
                        val: re.val,
                        node: node
                    }
                }
                else if (re.val_type == 'matrices') {
                    var height = 0;
                    var width = 0;
                    for (var i = 0; i < re.val.length; ++i) {
                        var mat = re.val[i]
                        height = mat.height > height ? mat.height : height
                        width = width + mat.width
                    }
                    var matrix_object = createMatrix('', height, width);
                    var bc = 0;
                    for (var i = 0; i < re.val.length; ++i) {
                        var mat = re.val[i];
                        var er = mat.height;
                        var ec = bc + mat.width;
                        for (var ir = 0; ir < er; ++ir) {
                            for (var ic = bc; ic < ec; ++ic) {
                                matrix_object.matrix[ir][ic] = mat.matrix[ir][ic - bc];
                            }
                        }
                        bc = ec;
                    }
                    return {
                        val_type: 'matrix',
                        val: matrix_object,
                        node: node
                    }
                }
            }
        }
        else if (node.mark[2] == '}') {
            return evalNode(node.l);
        }
        else if (node.mark[2] == ',') {
            var le = evalNode(node.l)
            var re = evalNode(node.r)
            if (le.val_type == 'number' && re.val_type == 'number') {
                return {
                    val: [le.val, re.val],
                    val_type: 'numbers',
                    node: node
                }
            }
            else if (le.val_type == 'number' && re.val_type == 'numbers') {
                return {
                    val: [le.val, ...re.val],
                    val_type: 'numbers',
                    node: node,
                }
            }
            else if (le.val_type == 'matrix' && re.val_type == 'matrix') {
                return {
                    val: [le.val, re.val],
                    val_type: 'matrices',
                    node: node,
                }
            }
            else if (le.val_type == 'matrix' && re.val_type == 'matrices') {
                return {
                    val: [le.val, ...re.val],
                    val_type: 'matrices',
                    node: node,
                }
            }
        }
        else if (node.mark[3] == 'name' && node.r == null) {
            var name = node.mark[2]
            var mat = findMatrixObjectByName(name)
            if (mat == null) {
                mat = createMatrix(name)
                matrices.push(mat)
            }
            return {
                val_type: 'matrix',
                val: mat,
                node: node
            }
        }
        else if (node.mark[3] == 'name' && node.r) {
            var re = evalNode(node.r);
            if (re.val_type == 'number') {
                var name = node.mark[2]
                var mat = findMatrixObjectByName(name);
                if (mat == null) {
                    mat = createMatrix(name)
                    matrices.push(mat)
                }
                return {
                    val_type: 'row',
                    val: mat.matrix[re.val],
                    node: node,
                }
            }
        }
        else if (node.mark[3] == 'number') {
            return {
                val_type: 'number',
                val: node.mark[2],
                node: node
            }
        }
        else if (node.mark[3] == 'method' && node.mark[2] == 'transpose') {
            var re = evalNode(node.r)
            if (re.val_type == 'matrix') {
                return {
                    val_type: 'matrix',
                    val: makeTranspose(re.val),
                    node: node
                }
            }
        }
        else if (node.mark[3] == 'method' && node.mark[2] == 'determinant') {
            var re = evalNode(node.r)
            if (re.val_type == 'matrix') {
                return {
                    val_type: 'number',
                    val: calcDeterminant(re.val),
                    node: node
                }
            }
        }
        else if (node.mark[3] == 'method' && node.mark[2] == 'minors') {
            var re = evalNode(node.r)
            if (re.val_type == 'matrix') {
                return {
                    val_type: 'matrix',
                    val: makeMinors(re.val),
                    node: node
                }
            }
        }
        else if (node.mark[3] == 'method' && node.mark[2] == 'cofactors') {
            var re = evalNode(node.r)
            if (re.val_type == 'matrix') {
                return {
                    val_type: 'matrix',
                    val: makeCofactors(re.val),
                    node: node
                }
            }
        }
        else if (node.mark[3] == 'method' && node.mark[2] == 'adjugate') {
            var re = evalNode(node.r)
            if (re.val_type == 'matrix') {
                return {
                    val_type: 'matrix',
                    val: makeTranspose(makeCofactors(re.val)),
                    node: node
                }
            }
        }
        else if (node.mark[3] == 'method' && node.mark[2] == 'inverse') {
            var re = evalNode(node.r)
            if (re.val_type == 'matrix') {
                return {
                    val_type: 'matrix',
                    val: makeInverse(re.val),
                    node: node
                }
            }
        }

        throw 'unhandled syntax error'
    }
    function displayResults(statements, results) {
        for (var i = 0; i < results.length; ++i) {
            var statement = statements[i]
            if (statement.length == 0)
                continue;
            var result = results[i]
            if (result == null)
                continue;
            var val = result.val
            var val_type = result.val_type
            if (val_type == 'empty')
                continue;
            if (result.node.mark[2] == '=' && val_type == 'matrix') {
                var widget = createDisplayMatrixWidget(val.name, val)
                element_matrix_display.appendChild(widget)
            }
            else if (val_type == 'matrix') {
                var name = statement.filter(mark => mark[3] != 'whitespace').map(mark => mark[2]).join('')
                var widget = createDisplayMatrixWidget(name, val)
                element_matrix_display.appendChild(widget)
            }
            else if (val_type == 'number') {
                var name = statement.filter(mark => mark[3] != 'whitespace').map(mark => mark[2]).join('')
                element_matrix_display.appendChild(createDisplayNumberWidget(name, val))
            }
        }
    }
    function createDisplayMatrixWidget(name, mat) {
        var cloned = template_matrix_display_entry.content.cloneNode(true)
        var div = document.createElement('div')
        div.appendChild(cloned)
        var element_matrix_display_entry_name = div.getElementsByClassName('matrix_display_entry_name')[0]
        var element_matrix_output = div.getElementsByClassName('matrix_output')[0]
        div.getElementsByClassName('matrix_display_entry_name')
        var table = fpCreateTable(mat.height, mat.width, (ir, ic, i) => {
            var label = document.createElement('label');
            label.innerHTML = mat.matrix[ir][ic];
            return label;
        })
        element_matrix_output.classList.add('update_table_on_recalc')
        element_matrix_output.setAttribute('matrix_name', name)
        element_matrix_display_entry_name.innerText = name
        element_matrix_output.appendChild(table)
        div.children[0].targetMatrix = mat;
        return div.children[0]
    }
    function createDisplayNumberWidget(name, number) {
        var cloned = template_number_display_entry.content.cloneNode(true)
        var div = document.createElement('div')
        div.appendChild(cloned)
        var element_matrix_display_entry_name = div.getElementsByClassName('matrix_display_entry_name')[0]
        var element_number_output = div.getElementsByClassName('number_output')[0]
        element_matrix_display_entry_name.innerText = name;
        element_number_output.innerText = number;
        return div.children[0]
    }
    function serializeState() {
        return {
            matrices: JSON.stringify(matrices.filter(matrix => matrix.input_node != null)),
            expressions: element_writer_section_textarea.value
        }
    }
    window.serializeState = serializeState;
    // state : object returend by serializeState
    function loadFromState(state) {
        if (state.matrices) {
            matrices.push(...JSON.parse(state.matrices))
            for (var i = 0; i < matrices.length; ++i) {
                element_matrix_area.insertBefore(createMatrixInputWidget(matrices[i]), element_add_widget_button)
            }
        }
        if (state.expressions) {
            element_writer_section_textarea.value = state.expressions;
        }
    }
    function saveStateToLocalStorage() {
        var state = serializeState();
        window.localStorage.setItem('matrices', state.matrices);
        window.localStorage.setItem('expressions', state.expressions);
    }
    function loadStateFromLocalStorage() {
        return {
            matrices: window.localStorage.getItem('matrices'),
            expressions: window.localStorage.getItem('expressions')
        }
    }
}

window.addEventListener('load', onLoad)
