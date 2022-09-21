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
function precision(a) 
{
    if (!isFinite(a)) return 0;
    var e = 1, p = 0;
    while (Math.round(a * e) / e !== a) { e *= 10; p++; }
    return p;
}
function copyMatrix(dst, src)
{
    dst.width = src.width
    dst.height = src.height
    dst.matrix = JSON.parse(JSON.stringify(src.matrix))
}
function subMatrix(mat, br, bc, er, ec)
{
    var matrix_object = createMatrix(mat.name);
    resizeMatrix(matrix_object, ec - bc, er - br);
    for( var ir = br; ir < er; ++ir )
    {
        for( var ic = bc; ic < ec; ++ic)
        {
            matrix_object.matrix[ir - br][ic - bc] = mat.matrix[ir][ic];
        }
    }
    return matrix_object;
}
function createMatrix(matrix_name) {
    return {
        name: matrix_name,
        type: 'matrix',
        width: 2,
        height: 2,
        matrix: [[1, 0], [0, 1]],
    };
}
function addMatrices(a, b) {
    var r = a.height > b.height ? a.height : b.height;
    var c = a.width > b.width ? a.width : b.width;
    var matrix_object = createMatrix('?')
    reserveMatrix(a, r, c);
    reserveMatrix(b, r, c);
    resizeMatrix(matrix_object, r, c)
    for (var ir = 0; ir < r; ++ir)
    {
        for (var ic = 0; ic < c; ++ic)
        {
            matrix_object.matrix[ir][ic] = a.matrix[ir][ic] + b.matrix[ir][ic]
        }
    }
    return matrix_object
}
function multMatrices(a, b) {
    if (a.width == b.height) {
        var matrix_object = createMatrix('');
        resizeMatrix(matrix_object, a.height, b.width);
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
        return null
    }
}
function scaleMatrix(k, mat)
{
    var matrix_object = createMatrix('')
    resizeMatrix(matrix_object, mat.height, mat.width)
    for (var ir = 0; ir < mat.height; ++ir) {
        for (var ic = 0; ic < mat.width; ++ic) {
            matrix_object.matrix[ir][ic] = k * mat.matrix[ir][ic]
        }
    }
    return matrix_object
}
function makeMinor(mat, r, c)
{
    var matrix_object = createMatrix('')
    resizeMatrix(matrix_object, mat.height-1, mat.width-1);
    var idstr = 0;
    for(var ir = 0; ir < mat.height; ++ir )
    {
        if(ir == r)
        {
            continue;
        }
        var idstc = 0;
        for( var ic = 0; ic < mat.width; ++ic)
        {
            if(c == ic)
                continue
            matrix_object.matrix[idstr][idstc++] = mat.matrix[ir][ic];
        }
        ++idstr;
    }
    return matrix_object;
}
function makeMinors(mat)
{
    if( mat.width != mat.height )
        throw 'makeMinors() : non-square matrices do not have minors'
    var matrix_object = createMatrix('');
    resizeMatrix(matrix_object, mat.height, mat.width);
    for( var ir= 0; ir < mat.height; ++ir )
    {
        for( var ic = 0; ic < mat.width; ++ic)
        {
            var minor = makeMinor(mat, ir, ic)
            matrix_object.matrix[ir][ic] = calcDeterminant(minor)
        }
    }
    return matrix_object;
}
function makeCofactors(mat)
{
    if( mat.width != mat.height )
        throw 'makeCofactors() : non-square matrices do not have cofactors'
    var matrix_object = createMatrix('');
    resizeMatrix(matrix_object, mat.height, mat.width);
    for( var ir= 0; ir < mat.height; ++ir )
    {
        var sign = ir % 2 == 0 ? 1 : -1;
        for( var ic = 0; ic < mat.width; ++ic)
        {
            var minor = makeMinor(mat, ir, ic)
            matrix_object.matrix[ir][ic] = sign * calcDeterminant(minor)
            sign *= -1;
        }
    }
    return matrix_object;
}
function makeInverse(mat)
{
    if( mat.width != mat.height)
        throw 'makeInverse() : non-square matrices do not have inverses'
    var adjugate = makeTranspose(makeCofactors(mat))
    var determinant = calcDeterminant(mat);
    return scaleMatrix(1 / determinant, adjugate);
}
function calcDeterminant(mat)
{
    if( mat.width != mat.height )
    {
        throw `calcDeterminant() : non-square matrices do not have determinants. (${mat.height}×${mat.width})`
    }
    if (mat.width == 1 && mat.height == 1)
    {
        return mat.matrix[0][0]
    }
    var sign = 1;
    var determinant = 0;
    for ( var ic = 0; ic < mat.width; ++ic )
    {
        var a = mat.matrix[0][ic]
        var minor = makeMinor(mat, 0, ic)
        var minor_determinant = calcDeterminant(minor)
        var c = sign * a * minor_determinant
        determinant += c
        sign *= -1
    }
    return determinant
}
function makeTranspose(mat)
{
    var matrix_object = createMatrix('')
    resizeMatrix(matrix_object, mat.width, mat.height)
    for( var ir = 0; ir < mat.height; ++ir)
    {
        for( var ic = 0; ic < mat.width; ++ic )
        {
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
function isWhitespaceChar(c)
{
    if (c == undefined) return false
    return /\s/.test(c)
}
function isOperatorChar(c)
{
    if (c == undefined) return false
    if (c == '=' || c == '+' || c == '*' || c == ';' || c =='(' || c== ')') return true
    return false
}
function isNameChar(c)
{
    if (c == undefined) return false
    return !isWhitespaceChar(c) && !isOperatorChar(c)
}
function markWhitespace(str, from)
{
    var i = from
    while(isWhitespaceChar(str[i++]));
    return [from, i-1, str.substr(from, i-from-1), 'whitespace']
}
function markName(str, from)
{
    var i = from
    while(isNameChar(str[i++]));
    return [from, i-1, str.substr(from, i-from-1), 'name']
}
function markMethod(str, from)
{
    for(var i = 0; i < method_names.length; ++i )
    {
        var method_name = method_names[i]
        var substr = str.substr(from, method_name.length)
        if( substr == method_name)
        {
            return [from, from + method_name.length, method_name, 'method'];
        }
    }
    return [from, from, null, 'method'];
}
function markOperator(str, from)
{
    if (isOperatorChar(str[from]))
        return [from, from+1, str[from], 'operator']
    else
        return [from, from, null, 'operator']
}
function markNumber(str, from)
{
    var float_regex = /[\+\-]?((\d+(\.\d*)?)|(\.\d+))([eE][\+\-]?\d*)?/g
    float_regex.lastIndex = from
    var ary = float_regex.exec(str);
    if ( ary == null || ary.index != from )
    {
        return [from, from, NaN, 'number']
    }
    var number = parseFloat(ary[0])
    if ( isNaN(number) )
    {
        console.warn('markNumber() : parseFloat() failed, even after regex match')
        return [from, from, NaN, 'number']
    }
    return [from, from + ary[0].length, number, 'number']
}
function markSomething(str, from)
{
    var mark;
    mark = markWhitespace(str, from)
    if( mark[0] != mark[1] )
        return mark
    mark = markNumber(str, from)
    if( mark[0] != mark[1])
        return mark
    mark = markMethod(str, from)
    if( mark[0] != mark[1] )
        return mark
    mark = markOperator(str, from)
    if( mark[0] != mark[1] )
        return mark
    mark = markName(str, from)
    if( mark[0] != mark[1] )
        return mark
    return null
}
function markAll(str)
{
    if( str == null )
        return null;
    i = 0
    marks = []
    while( i != str.length )
    {
        mark = markSomething(str, i)
        if(mark[0] == mark[1])
            return null; // fail
        marks.push(mark)
        i = mark[1]
    }
    return marks
}
function onLoad(event)
{
    var template_matrix_input = document.getElementById("template_matrix_input")
    var template_matrix_display_entry = document.getElementById("template_matrix_display_entry")
    var element_matrix_area = document.getElementById("matrix_area")
    var element_add_widget_button = document.getElementById("add_widget_button")
    var element_writer_section_textarea = document.getElementsByClassName('writer_section_textarea')[0]
    var element_matrix_display = document.getElementById('matrix_display')
    var element_menu_exporter = document.getElementById('menu_exporter')
    element_writer_section_textarea.addEventListener('keyup', (event) =>
    {
        if( event.target == element_writer_section_textarea )
        {
            recalcMatrix()
        }
    })
    loadState();
    function evalStatement(marks)
    {
        nodes = marks.filter(mark => mark[3] != 'whitespace').map(mark => {
            return {
                l: null,
                r: null,
                mark: mark,
                val: null,
            }
        })
        if( nodes.length > 0 )
        {
            root = treeficate(nodes)
            var result = evalNode(root)
            
        }
    }
    function renewMatrixList(marks)
    {
        var marks = marks.filter(mark=>mark[3] == 'name')
        for( var i = 0; i < marks.length; ++i )
        {
            var mark = marks[i];
            var name = mark[2];
            var mo = matrices.findIndex(m => m.name == name)
            if ( mo == -1 )
            {
                matrices.push(createMatrix(name))
            }
        }
        for( var i = matrices.length - 1; i >= 0; --i)
        {
            var mo = matrices[i]
            if (marks.findIndex(mark=>mark.name == mo.name) == -1 && !mo.input_node)
            {
                matrices.splice(i, 1);
            }
        }
    }
    function evalScript()
    {
        str = element_writer_section_textarea.value // something like '   C = A + B   ...'
        if( str == '' )
            return;
        var marks = markAll(str)
        renewMatrixList(marks)
        var statements = divideStatements(marks)
        for( var i = 0; i < statements.length; ++i )
        {
            var statement = statements[i]
            evalStatement(statement)
        }
        displayMatrices()
    }
    function recalcMatrix()
    {
        evalScript()
        var updates = document.getElementsByClassName('update_table_on_recalc')
        for (var i = 0; i < updates.length; ++i )
        {
            var name = updates[i].getAttribute("matrix_name")
            var mat = findMatrixObjectByName(name)
            if( mat == null )
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
    }
    function createMatrixInputNode(ir, ic, i, matrix_object)
    {
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
    function establishMatrixInput(matrix_input_node, matrix_object)
    {
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
    function createMatrixInputWidget(matrix_object)
    {
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
    window.onclick = function (event)
    {
        if (event.target.matches('.export_button > svg'))
        {
            event.target.parentElement.appendChild(element_menu_exporter)
            element_menu_exporter.classList.add('show')
        }
        else
        {
            if (element_menu_exporter.classList.contains('show')) {
                element_menu_exporter.classList.remove('show');
            }
        }
        if( event.target.matches('#a_json_row_major'))
        {
            var widget = event.target.parentElement.parentElement.parentElement;
            var element_matrix_output = widget.getElementsByClassName('matrix_output')[0]
            var name = element_matrix_output.getAttribute("matrix_name")
            var matrix = matrices.find(m=>m.name == name)
            var to_stringify = subMatrix(matrix, 0, 0, matrix.height ,matrix.width);
            navigator.clipboard.writeText(JSON.stringify(to_stringify.matrix))
            event.preventDefault()
        }
        if( event.target.matches('#a_json_col_major'))
        {
            var widget = element_menu_exporter.parentElement.parentElement;
            var element_matrix_output = widget.getElementsByClassName('matrix_output')[0]
            var name = element_matrix_output.getAttribute("matrix_name")
            var matrix = matrices.find(m=>m.name == name)
            navigator.clipboard.writeText(JSON.stringify(makeTranspose(subMatrix(matrix, 0, 0, matrix.height, matrix.width)).matrix))
            event.preventDefault()
        }
        if( event.target.matches('#a_katex'))
        {
            var widget = element_menu_exporter.parentElement.parentElement;
            var element_matrix_output = widget.getElementsByClassName('matrix_output')[0]
            var name = element_matrix_output.getAttribute("matrix_name")
            var matrix = matrices.find(m=>m.name == name)
            var str = '\\begin{pmatrix}'
            for( var ir = 0; ir < matrix.height; ++ir )
            {
                for( var ic = 0; ic < matrix.width; ++ic )
                {
                    var element = matrix.matrix[ir][ic]
                    str += element.toString()
                    if( ic < matrix.width - 1 )
                        str += '&'
                    else
                        str += '\\\\'
                }
            }
            str += '\\end{pmatrix}'
            navigator.clipboard.writeText(str)
            event.preventDefault()
        }
        if( event.target.matches('#a_c_row_major'))
        {
            var widget = element_menu_exporter.parentElement.parentElement;
            var element_matrix_output = widget.getElementsByClassName('matrix_output')[0]
            var name = element_matrix_output.getAttribute("matrix_name")
            var matrix = matrices.find(m=>m.name == name)
            var str = JSON.stringify(matrix.matrix).replaceAll('[','{').replaceAll(']','}')
            navigator.clipboard.writeText(str)
            event.preventDefault()
        }
        if( event.target.matches('#a_c_col_major'))
        {
            var widget = element_menu_exporter.parentElement.parentElement;
            var element_matrix_output = widget.getElementsByClassName('matrix_output')[0]
            var name = element_matrix_output.getAttribute("matrix_name")
            var matrix = makeTranspose(matrices.find(m=>m.name == name))
            var str = JSON.stringify(matrix.matrix).replaceAll('[','{').replaceAll(']','}')
            navigator.clipboard.writeText(str)
            event.preventDefault()
        }
        if(event.target.matches('#add_widget_button'))
        {
            var names = matrices.map(matrix_object => matrix_object.name)
            
            function nextChar(char)
            {
                if( char == 'H') return 'J';
                if( char == 'N') return 'P';
                return String.fromCharCode(char.charCodeAt(0) + 1);
            }
            function nextName(name)
            {
                var new_name_ary = []
                var carry = 1;
                for( var i = name.length - 1; i >= 0; --i)
                {
                    if( carry == 0)
                    {
                        new_name_ary.push(name[i])
                    }
                    else
                    {
                        if( name[i] == 'Z')
                        {
                            carry = 1
                            new_name_ary.push('A')
                        }
                        else
                        {
                            carry = 0
                            new_name_ary.push(nextChar(name[i]));
                        }
                    }
                }
                if( carry == 1 )
                {
                    new_name_ary.push('A')
                }
                return new_name_ary.reverse().join('')
            }
            var new_name = 'A';
            while( names.indexOf(new_name) != -1)
            {
                new_name = nextName(new_name);
            }
            var matrix_object = createMatrix(new_name)
            matrices.push(matrix_object)
            element_matrix_area.insertBefore(createMatrixInputWidget(matrix_object), element_add_widget_button)
            event.preventDefault()
        }
    }
    function divideStatements(marks)
    {
        var statements = []
        var statement = []
        for( var i = 0; i < marks.length; ++i)
        {
            var mark = marks[i]
            if( mark[2] == ';' )
            {
                statements.push(statement)
                statement = []
                continue;
            }
            else
            {
                statement.push(marks[i])
            }
        }
        return statements
    }
    function treeficate(nodes)
    {
        var to_fold = []
        var depth = 0;
        for( var i = 0; i < nodes.length; ++i )
        {
            var node = nodes[i]
            if(node.mark[2] == '(')
            {
                nodes.splice(i, 1)
                --i;
                ++depth;
                continue;
            }
            if(node.mark[2] == ')')
            {
                nodes.splice(i, 1)
                --i;
                --depth;
                continue;
            }
            if(node.mark[3] == 'method')
            {
                to_fold.push({
                    node : node,
                    depth : depth
                })
            }
            if(node.mark[3] == 'operator')
            {
                to_fold.push({
                    node : node,
                    depth : depth
                })
            }
        }
        var optier =
        {
            '=':0,
            '+':1,
            '-':1,
            '*':2,
            '/':2,
            'method':3,
        };
        fpHash = function (a)
        {
            var str = a.node.mark[2]
            var type = a.node.mark[3]
            if( type == 'operator' )
            {
                return a.depth * 100 + optier[str]
            }
            else if (type == 'method')
            {
                return a.depth * 100 + optier['method']
            }
            throw 'type must be operator|method'
        }
        to_fold.sort((a, b) => {
            return fpHash(a) - fpHash(b)
        })
        for( var d = to_fold.length-1; d >= 0; --d)
        {
            var fold = to_fold[d]
            var opnode = fold.node
            var index = nodes.indexOf(opnode)
            if (opnode.mark[3] == 'operator')
            {
                opnode.r = nodes.splice(index+1, 1)[0];
                opnode.l = nodes.splice(index-1, 1)[0];
            }
            if (opnode.mark[3] == 'method')
            {
                opnode.r = nodes.splice(index+1, 1)[0];
            }
        }
        return nodes[0]
    }
    function evalNode(node)
    {
        if( node == null )
            throw 'evalNode() : node is null'
        if( node.mark == null )
            throw 'evalNode() : node.mark is null'

        if( node.mark[2] == '=' )
        {
            var le = evalNode(node.l)
            var re = evalNode(node.r)
            if( le.val_type == 'matrix' && re.val_type =='matrix' )
            {
                copyMatrix(le.val, re.val)
                return {
                    val_type : 'matrix',
                    val: le,
                }
            }
        }
        else if( node.mark[2] == '+')
        {
            var le = evalNode(node.l)
            var re = evalNode(node.r)
            if( le.val_type == 'matrix' && re.val_type == 'matrix')
            {
                return {
                    val_type : 'matrix',
                    val : addMatrices(le.val, re.val)
                }
            }   
            else if( le.val_type == 'number' && re.val_type == 'number')
            {
                return {
                    val_type : 'number',
                    val : le.val + re.val
                }
            }
        }
        else if( node.mark[2] == '*')
        {
            var le = evalNode(node.l)
            var re = evalNode(node.r)
            if( le.val_type == 'matrix' && re.val_type == 'matrix')
            {
                return {
                    val_type : 'matrix',
                    val : multMatrices(le.val, re.val)
                }
            }
            else if( le.val_type == 'number' && re.val_type == 'matrix')
            {
                return {
                    val_type : 'matrix',
                    val : scaleMatrix(le.val, re.val)
                }
            }
            else if( le.val_type == 'matrix' && re.val_type == 'number')
            {
                return {
                    val_type : 'matrix',
                    val : scaleMatrix(re.val, le.val)
                }
            }
            else if( le.val_type == 'number' && re.val_type == 'number')
            {
                return {
                    val_type : 'number',
                    val : le.val * re.val
                }
            }
        }
        else if( node.mark[3] == 'name')
        {
            var name = node.mark[2]
            var mat = findMatrixObjectByName(name)
            if (mat == null )
            {
                mat = createMatrix(name)
                matrices.push(mat)
            }
            return {
                val_type : 'matrix',
                val : mat
            }
        }
        else if( node.mark[3] == 'number' )
        {
            return {
                val_type : 'number',
                val : node.mark[2]
            }
        }
        else if( node.mark[3] == 'method' && node.mark[2] == 'transpose' )
        {
            var re = evalNode(node.r)
            if( re.val_type == 'matrix' )
            {
                return {
                    val_type : 'matrix',
                    val : makeTranspose(re.val)
                }
            }
        }
        else if( node.mark[3] == 'method' && node.mark[2] == 'determinant')
        {
            var re = evalNode(node.r)
            if( re.val_type == 'matrix')
            {
                return{
                    val_type : 'number',
                    val : calcDeterminant(re.val)
                }
            }
        }
        else if(node.mark[3] == 'method' && node.mark[2] == 'minors')
        {
            var re = evalNode(node.r)
            if( re.val_type == 'matrix')
            {
                return {
                    val_type : 'matrix',
                    val : makeMinors(re.val)
                }
            }
        }
        else if(node.mark[3] == 'method' && node.mark[2] == 'cofactors')
        {
            var re = evalNode(node.r)
            if( re.val_type == 'matrix')
            {
                return {
                    val_type : 'matrix',
                    val : makeCofactors(re.val)
                }
            }
        }
        else if(node.mark[3] == 'method' && node.mark[2] == 'adjugate')
        {
            var re = evalNode(node.r)
            if( re.val_type == 'matrix')
            {
                return {
                    val_type : 'matrix',
                    val : makeTranspose(makeCofactors(re.val))
                }
            }
        }
        else if(node.mark[3] == 'method' && node.mark[2] == 'inverse')
        {
            var re = evalNode(node.r)
            if( re.val_type == 'matrix')
            {
                return {
                    val_type : 'matrix',
                    val : makeInverse(re.val)
                }
            }
        }
    }
    function displayMatrices()
    {
        element_matrix_display.innerText = ''
        for( var i = 0; i < matrices.length; ++i )
        {
            if( matrices[i].input_node != null)
                continue;
            var cloned = template_matrix_display_entry.content.cloneNode(true)
            var div = document.createElement('div')
            div.appendChild(cloned)
            var mat = matrices[i]
            var element_matrix_display_entry_name = div.getElementsByClassName('matrix_display_entry_name')[0]
            var element_matrix_output = div.getElementsByClassName('matrix_output')[0]
            div.getElementsByClassName('matrix_display_entry_name')
            var table = fpCreateTable(mat.height, mat.width, (ir, ic, i) => {
                var label = document.createElement('label');
                label.innerHTML = mat.matrix[ir][ic];
                return label;
            })
            element_matrix_output.classList.add('update_table_on_recalc')
            element_matrix_output.setAttribute('matrix_name',mat.name)
            element_matrix_display_entry_name.innerText = mat.name
            element_matrix_output.appendChild(table)
            element_matrix_display.appendChild(div.children[0])
        }
    }
    function saveState()
    {
        var str = JSON.stringify(matrices.filter(matrix => matrix.input_node != null));
        window.localStorage.setItem('matrices', str)

        var str = element_writer_section_textarea.value
        window.localStorage.setItem('expressions', str)
    }
    function loadState()
    {
        var ls_matrices = window.localStorage.getItem('matrices')
        if( ls_matrices )
        {
            matrices.push(...JSON.parse(ls_matrices))
            for( var i = 0; i < matrices.length; ++i )
            {
                element_matrix_area.insertBefore(createMatrixInputWidget(matrices[i]), element_add_widget_button)
            }
        }
        var ls_expressions =window.localStorage.getItem('expressions')
        element_writer_section_textarea.value = ls_expressions;
    }
    function onBeforeUnload(event)
    {
        saveState()
    }
    window.addEventListener('beforeunload', onBeforeUnload)
}

window.addEventListener('load', onLoad)
