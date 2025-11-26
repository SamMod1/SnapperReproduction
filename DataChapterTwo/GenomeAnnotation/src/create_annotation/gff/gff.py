from pandas import DataFrame
from .features import Features
from create_annotation.general_funcs import progress_bar


def build_feature_fields(split_line):
    feature = {
        'sequence': split_line[0],
        'source': split_line[1],
        'feature_type': split_line[2],
        'start': int(split_line[3]),
        'end': int(split_line[4]),
        'score': split_line[5],
        'strand': split_line[6],
        'phase': split_line[7],
    }
    
    for attribute in split_line[8].split(';'):
        attribute_name, attribute = attribute.split('=', 1)
        feature[attribute_name] = attribute.replace('\n', '')
    
    return feature
                        

class Gff:
    def __init__(self, filename = None, filelines = None, data = None):
        self.data = None
        
        if not data is None:
            self.data = data
            
        elif filename is None and filelines is None and data is None:
            self.data = Features()
            
        else:
            if filename:
                origional_gff = open(filename, 'r').readlines()
            else:
                origional_gff = filelines
                
            self.data = Features()
            self.data.append({'ID': '_comments_', '_children_': [], 'feature_type': 'comments index'})
            
            comment_count = 0
            for line in origional_gff:
                if not line or line[0] == '#':
                    self.data.append({'ID': f'_comment{str(comment_count)}_', 'Parent': '_comments_', 'feature_type': 'comment', 'comment': line.replace('\n', '')})
                    comment_count += 1
                    continue
                line = line.split('\t')
                    
                feature = build_feature_fields(line)
                self.append(feature)

                
    def __getitem__(self, key):
        key_data = self.data[key]
        new_data = Features()

        if 'ID' in key_data:
            key_data = {key_data['ID']: key_data}
            
        for x in key_data:
            # if 'Parent' in key_data[x]:
            #     new_data.append(self.data[key_data[x]['Parent']])
            new_data.append(key_data[x])
            if '_children_' in key_data[x]:
                for child in key_data[x]['_children_']:
                    new_data.append(self.data[child])
                    
        return Gff(data = new_data)
                    
    def __iter__(self):
        return (x for x in self.data)
        
    def __str__(self):
        if len(self) < 10:
            return str(self.data.features)
        else:
            to_return = 'Feature Counts:\n\n'
            feature_counts = self.feature_counts()
            for x in feature_counts:
                to_return += x + ': ' + str(feature_counts[x]) + '\n'
                
            return to_return[:-1]
    
    def __len__(self):
        try:
            return len(self.data) - (len(self.data['_comments_']['_children_']) + 1)
        except KeyError:
            return len(self.data)
        
    def __setitem__(self, position, item):
        pass
    
    def get_feature(self, key, del_children = True):
        to_return = self.data.return_individual(key, del_children)
        return to_return
    
    def get_entire_lineage(self, key):
        while True:
            feature = self.get_feature(key)
            if 'Parent' in feature:
                return self.get_entire_lineage(feature['Parent'])
            else:
                return self[key]
    
    def find_replace(self, feature_name, replace, find = None):
        if not find:
            find = feature_name
        lineage = self.get_entire_lineage(feature_name)
        for feature_name in lineage:
            feature = lineage.get_feature(feature_name, del_children = False)
            for attribute_name in feature:
                if attribute_name == '_children_':
                    new_children = []
                    for child in feature['_children_']:
                        new_children.append(child.replace(find, replace))
                    self.edit_feature_data(feature['ID'], '_children_', new_children)
                elif type(feature[attribute_name]) == str:
                    if find in feature[attribute_name]:
                        if attribute_name == 'ID':
                            new_id = feature['ID'].replace(find, replace)
                            self.data.rename_single_feature(feature['ID'], new_id)
                            feature['ID'] = new_id
                        else:
                            self.edit_feature_data(feature['ID'], attribute_name, feature[attribute_name].replace(find, replace))

        #self.data.rename_single_feature(old_name, new_name)
                        
    def feature_counts(self):
        feature_counts = {'sequences': 0}
        sequences = []
        for x in self.data:
            x = self.data[x]
            try:
                if 'feature_type' in x:
                    feature_counts[x['feature_type']] += 1
                    if not x['sequence'] in sequences:
                        feature_counts['sequences'] += 1
                        sequences.append(x['sequence'])
            except KeyError:
                feature_counts[x['feature_type']] = 1
                    
        return feature_counts
    
    def append(self, feature):
        self.data.append(feature)
        
    def extend(self, second_gff):
        self.insert(second_gff, len(self.data))
        
    def insert(self, features, idx):
        if type(features) == dict:
            self.data.insert(features, idx)
        else:
            for count, feature in enumerate(features):
                self.insert(features.data[feature], idx + count)
                
    def pop(self, feature):
        features_to_remove = self[feature]
        if type(features_to_remove) == dict:
            return self.data.pop(feature)
        
        to_return = Gff(data=Features())
        for feature in features_to_remove:
            to_return.append(self.data.pop(feature))
        
        return to_return
                
    def to_df(self):
        try:
            comments = self.data['_comments_']['_children_']
            if not comments:
                comments = ['_comments_']
        except KeyError:
            comments = ['_comments_']
            
        df = {}
        for idx, feature in enumerate(self.data):
            if not feature in comments:
                feature = dict(self.data[feature])
                for existing_key in df:
                    if not existing_key in feature.keys():
                        feature[existing_key] = None
                
                for key in feature:
                    try:
                        df[key].append(feature[key])
                    except KeyError:
                        df[key] = [None for _ in range(idx)]
                        df[key].append(feature[key])
               
        if '_children_' in df:
            del df['_children_']
                        
        df = DataFrame(df)
        
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)
        
        return df
    
    def to_gff(self, show_progress_bar = False):        
        try:
            comments = self.data['_comments_']['_children_'] + ['_comments_']
            if not comments:
                comments = ['_comments_']
        except KeyError:
            comments = ['_comments_']
        
        self_len = len(self)
        non_attributes = set([
            'sequence', 'source', 'feature_type', 'start', 'end', 
            'score', 'strand', 'phase', '_children_', '_count_'
        ])
        lines = []
        line = ''
        for idx, feature_name in enumerate(self.data):
            
            if show_progress_bar:
                print('\r' + progress_bar(idx, self_len), end='\r')
            
            feature = dict(self.data[feature_name])
            if feature_name in comments:
                if not feature_name == '_comments_':
                    line += feature['comment'] + '\n'
            else:
                line += feature['sequence'] + '\t'
                line += feature['source'] + '\t'
                line += feature['feature_type'] + '\t'
                line += str(feature['start']) + '\t'
                line += str(feature['end']) + '\t'
                line += feature['score'] + '\t'
                line += feature['strand'] + '\t'
                line += feature['phase'] + '\t'
                
                for attribute in feature:
                    if not attribute in non_attributes:
                        line += attribute + '=' + str(feature[attribute]) + ';'
                    
                line = line[:-1]
                line += '\n'
                
            lines.append(line)
            line = ''
                
        return ''.join(lines)[:-1]
    
    def index(self, item):
        return self.data.index(item)
    
    def edit_feature_data(self, feature, to_edit, new_value):
        self.data[feature][to_edit] = new_value
        
    def add_attribute(self, feature, attribute_name, attribute_value):
        self.data[feature][attribute_name] = attribute_value
        
    def copy(self):
        new_data = self.data.copy()
        return Gff(data = new_data)
    
    def remove_single_feature(self, feature):
        self.data.remove_single_feature(feature)
        
    def edit_attributes(self, attribute_name, new_data):
        if type(new_data) == str:
            new_data = [new_data for _ in range(len(self))]
            
        for idx, feature in enumerate(self):
            feature = self.data[feature]
            feature[attribute_name] = new_data[idx]
                        
    def strip_comments(self):
        self.pop('_comments_')
