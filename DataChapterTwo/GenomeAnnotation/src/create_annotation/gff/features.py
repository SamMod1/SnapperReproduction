from copy import deepcopy


class Features:
    def __init__(self, features = None, feature_index = None):
        if not features is None and not feature_index is None:
            self.features = features
            self.feature_index = feature_index
        else:
            self.features = {}
            self.feature_index = []
    
    def __getitem__(self, key):
        if type(key) == str:
            return self.features[key]
        elif type(key) == int:
            return self.features[self.feature_index[key]]
        elif type(key) == slice:
            key = self.feature_index[key]
            
        to_return = Features()

        for item in key:
            to_return.append(self.features[item])

        return to_return
    
    def __len__(self):
        return len(self.feature_index)
    
    def __iter__(self):
        return (x for x in self.feature_index)
    
    def __setitem__(self, item):
        exists = item['ID'] in self.feature_index
        
        if not exists:
            self.append(item)
            if 'Parent' in item:
                self.add_child(item['ID'], item['Parent'])
        else:
            if 'Parent' in item:
                self.add_child(item['ID'], item['Parent'])

            idx = self.feature_index.index(item['ID'])
            self.features[item['ID']] = item
            self.feature_index.insert(idx, item['ID'])

    def append(self, item):
        if item['ID'] in self.features:
            try:
                self.features[item['ID']]['_count_'] += 1
                item['ID'] = item['ID'] + '-' + str(self.features[item['ID']]['_count_'])
            except KeyError:
                self.features[item['ID']]['_count_'] = 2
                item['ID'] = item['ID'] + '-2'
        self.features[item['ID']] = item
        self.feature_index.append(item['ID'])
        if 'Parent' in item:
            self.add_child(item['ID'], item['Parent'])
        
    def insert(self, item, idx):
        if type(item) == dict:
            self.features[item['ID']] = item
            self.feature_index.insert(idx, item['ID'])
            if 'Parent' in item:
                self.add_child(item['ID'], item['Parent'])
        else:
            for count, feature in enumerate(item):
                self.insert(feature, idx + count)
                
    def pop(self, to_pop):
        feature = dict(self.features[to_pop])
        if 'Parent' in feature:
            self.remove_child(feature['ID'], feature['Parent'])
        del self.features[to_pop]
        self.feature_index.remove(to_pop)
        
        return feature
        
    def add_child(self, child, parent):
        try:
            if '_children_' in self.features[parent]:
                if not child in self.features[parent]['_children_']:
                    self.features[parent]['_children_'].append(child)
            else:
                self.features[parent]['_children_'] = [child]

        except KeyError:
            return#self.features[parent]['_children_'] = [child]
            
        if 'Parent' in self.features[parent]:
            self.add_child(child, self.features[parent]['Parent'])
            
    def remove_child(self, child, parent):
        try:
            self.features[parent]['_children_'].remove(child)
                
            if 'Parent' in self.features[parent]:
                self.remove_child(child, self.features[parent]['Parent'])
        except KeyError:
            pass
        
    def index(self, item):
        return self.feature_index.index(item)
    
    def return_individual(self, key, del_children = True):
        to_return = deepcopy(self.features[key])
        if '_children_' in to_return and del_children:
            del to_return['_children_']
        return to_return
    
    def copy(self):
        new_idx = list(self.feature_index)
        new_features = {}
        for key in self.features:
            new_features[key] = dict(self.features[key])
            
        return Features(features = new_features, feature_index = new_idx)
    
    def remove_single_feature(self, feature):
        del self.features[feature]
        self.feature_index.remove(feature)
        
    def rename_single_feature(self, old_name, new_name, count = 1):
        if new_name in self.features:
            raise KeyError(f'{new_name} is already the name of a feature')
        feature_data = self.features[old_name]
        feature_index = self.feature_index.index(old_name)
        feature_data['ID'] = new_name
        self.remove_single_feature(old_name)
        self.insert(feature_data, feature_index)
        
        