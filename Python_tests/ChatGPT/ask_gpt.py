import pickle

file = open("stored_data.pickle", "rb")
all_data = pickle.load(file)
file.close()

docsearch = all_data["docsearch"]
chain = all_data["chain"]

requete = ""

while requete != "exit":
    requete = input("> ")
    documents_similaires = docsearch.similarity_search(requete)
    reponse = chain.run(input_documents=documents_similaires, question=requete)
    print(reponse)
